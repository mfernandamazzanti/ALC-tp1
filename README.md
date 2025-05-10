# TP1 – Museum Network in Buenos Aires  
## Computational Linear Algebra – UBA 2025

This repository contains the first group assignment for the course **Computational Linear Algebra** at the University of Buenos Aires (FCEN, UBA – 1st semester, 2025).

---

## 📌 Overview

The goal of this project is to model a network of museums in Buenos Aires using graph theory and apply linear algebra techniques to analyze its structure and behavior.

We work with:
- **Adjacency and transition matrices**
- **PageRank algorithm**
- **LU factorization (implemented from scratch)**
- **Condition number analysis**
- **Geographic visualizations** of the museum network

---

## 📁 Files

- `TP1.ipynb`: Jupyter notebook with theoretical development, Python code, visualizations, and written answers.
- `template_funciones.py`: Python module containing manually implemented functions (no direct matrix inversion).

---

## 🗺️ Visualizations

Using `geopandas`, `matplotlib`, and `networkx`, the project includes:

- A map of Buenos Aires with museum nodes  
- Node sizes scaled by PageRank values  
- Comparative visualizations for different values of:
  - `m`: number of neighbor connections per museum  
  - `α`: damping factor (restart probability)



