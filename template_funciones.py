# Carga de paquetes necesarios para graficar
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd # Para leer archivos
import geopandas as gpd # Para hacer cosas geográficas
import seaborn as sns # Para hacer plots lindos
import networkx as nx # Construcción de la red en NetworkX
import scipy

#Funciones auxiliares (usadas en ejemplos y en otras funciones)
def calcula_matriz_K(A):
    conexiones = A.sum(axis=1)  # suma de cada fila (axis=1), devuelve un array
    K = np.diag(conexiones) #np diag crea una matriz diagonal a partir de una lista de elementos
    return K

def calcula_matriz_C(A):
    # Función para calcular la matriz de trancisiones C
    # A: Matriz de adyacencia
    # Retorna la matriz C
    A = np.array(A) #por si A viene en series de pandas, la pasamos a arrays de numpy
    K = calcula_matriz_K(A)
    Atrasp = A.T # Transponemos la matriz A
    Kinv = np.diag(1/np.diag(K)) # Como K es una matriz diagonal, la inversa es simplemente la inversa de los elementos de la diagonal
    C = Atrasp @ Kinv
    return C

def calcula_matriz_M(A, alpha):
    A = np.array(A) #por si la matriz esta en listas de python la paso a arryas de numpy
    C = calcula_matriz_C(A) #calculo la matriz C de A con la funcion definida previamente
    N = A.shape[0] #calculo el tamañao de las filas de A
    I = np.eye(N) #genero una matriz identidad de ese tamaño
    M = (N/alpha) * (I - (1 - alpha )*C) #genero la matriz M segun su definicion
    return M #devuelve la matriz M de A


#Funciones principales
def construye_adyacencia(D,m): 
    # Función que construye la matriz de adyacencia del grafo de museos
    # D matriz de distancias, m cantidad de links por nodo
    # Retorna la matriz de adyacencia como un numpy.
    D = D.copy()
    l = [] # Lista para guardar las filas
    for fila in D: # recorriendo las filas, anexamos vectores lógicos
        l.append(fila<=fila[np.argsort(fila)[m]] ) # En realidad, elegimos todos los nodos que estén a una distancia menor o igual a la del m-esimo más cercano
    A = np.asarray(l).astype(int) # Convertimos a entero
    np.fill_diagonal(A,0) # Borramos diagonal para eliminar autolinks
    return(A)

#Ej 3
def calculaLU(matriz):
    # matriz es una matriz de NxN
    # Retorna la factorización LU a través de una lista con dos matrices L y U de NxN.
    m=matriz.shape[0]
    n=matriz.shape[1]
    Ac = matriz.astype(float).copy() # Convertimos a float para evitar problemas de precisión

    if m!=n:
        print('Matriz no cuadrada')
        return

    for j in range(m): # Recorremos las columnas
      assert(Ac[j,j]!=0) # Usamos este assert para evitar divisiones por cero
      for i in range(j+1,n): # Recorremos las filas (abajo de la diagonal)
        Ac[i,j] = Ac[i,j]/Ac[j,j] #Calculamos el multiplicador
        for k in range(j+1,n): # Recorremos las columnas (a la derecha de la diagonal) para obtener U
          Ac[i,k] = Ac[i,k] - Ac[i,j]*Ac[j,k] # Restamos el multiplicador por la fila correspondiente de U

    L = np.tril(Ac,-1) + np.eye(matriz.shape[0]) # Matriz L es la parte inferior de Ac + matriz identidad
    U = np.triu(Ac) # Matriz U es la parte superior de Ac

    return L, U

def calcula_pagerank(A,alfa):
    # Función para calcular PageRank usando LU
    # A: Matriz de adyacencia
    # alfa: parámetro de amortiguación (damping factor), 0 < alfa < 1
    # Retorna: Un vector p con los coeficientes de page rank de cada museo
    C = calcula_matriz_C(A)
    N = A.shape[0] #tamañao de la primera fila=cantidad de museos
    M = calcula_matriz_M(A,alfa) #calculamos la matriz M de A
    L, U = calculaLU(M) # Calculamos descomposicion LU a partir de C y alfa
    b = np.ones(N) # Vector de 1s, multiplicado por el coeficiente correspondiente usando d y N.
    Up = scipy.linalg.solve_triangular(L,b,lower=True) # Primera inversion usando L
    p = scipy.linalg.solve_triangular(U,Up) # Segunda inversion usando U
    return p



#Ej 5
def calcula_norma_1(v):
    # Función para calcular la norma 1 de un vector 
    norma = 0
    for elem in v: 
        norma += abs(elem)
    return norma

def calcula_matriz_C_continua(D): 
    # Función para calcular la matriz de trancisiones C
    # D: matriz de distancias entre museos 
    # Retorna la matriz C en versión continua
    D = D.copy()
    np.fill_diagonal(D,np.inf) # evita división por cero en la diagonal
    F = 1/D
    np.fill_diagonal(F,0)
    K = calcula_matriz_K(F) 
    Kinv = np.diag(1/np.diag(K)) # Como K es una matriz diagonal, la inversa es simplemente la inversa de los elementos de la diagonal
    C = F @ Kinv # Calcula C multiplicando Kinv y F
    return C

def calcula_B(C,r):
    # C: Matirz de transiciones
    # r: Cantidad de pasos en la red dado por los visitantes. Indicado como r en el enunciado
    # Retorna: Una matriz B que vincula la cantidad de visitas w con la cantidad de primeras visitas v
    B = np.eye(C.shape[0]) #Voy guardando en B el resultado de la sumatoria (ya está C^0 = I)
    Ck = np.eye(C.shape[0]) #Inicializo C^k en C^0 = I
    for i in range(1,r):
        Ck = Ck @ C
        B += Ck 
    return B


#Ej 6
def calcular_inversa(A): #La idea es resolver el sistema AX = I, donde I es la matriz identidad y X será la matriz inversa de A
    L, U = calculaLU(A)
    I = np.eye(A.shape[0])
    Y = scipy.linalg.solve_triangular(L, I, lower=True)  # Resuelve LY = I
    X = scipy.linalg.solve_triangular(U, Y, lower=False)  # Resuelve UX = Y
    return X  # Matriz inversa de A

def calcula_norma_inducida_1(A):
    # Función para calcular la norma inducida 1 de una matriz
    filas = A.shape[0] #cantidad de filas de la matriz
    columnas = A.shape[1] #cantidad de columnas de la matriz
    columna_maxima = 0
    for j in range(columnas):
        suma_col = 0
        for i in range(filas):
            suma_col += abs(A[i][j])
        if suma_col > columna_maxima:
            columna_maxima = suma_col
    return columna_maxima #devuelve la suma máxima de columnas

def calcula_condicion_1_B(B):
  # Función para calcular el número de condición 1 de la matriz B
  norma_B = calcula_norma_inducida_1(B) #calculo la norma 1 de B
  Binv = calcular_inversa(B) #calculo la inversa de B con la funcion que armamos para calcular inversas con la factorizacion LU
  norma_B_inv = calcula_norma_inducida_1(Binv) #calculo la norma 1 de B inversa
  condicion = norma_B * norma_B_inv
  return condicion

