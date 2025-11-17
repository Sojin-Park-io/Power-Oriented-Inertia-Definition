# Power-Oriented-Inertia-Definition

Simulation data and code for the paper titled: **"Definition and Formulation of Inertia Service Incorporating Inverter-Based Resources."**

This repository contains the simulation code and data used to define and formalize inertia service incorporating inverter-based resources as presented in the paper.

---

## 🛠️ Development Environment and Requirements

The following setup was used to perform the simulations and solve the optimization problems presented in this research.

### 💻 Hardware Specifications

All simulations were performed on a laptop computer with the following specifications:
* **CPU:** Apple M3 chip (4.05 GHz clock rate)
* **RAM:** 16 GB

### ⚙️ Software and Tools

To run the code and solve the optimization problems, the following software and tools must be installed:

* **Main Environment:** **Matlab (R2025a)**
* **Base Dataset Source:** **MATPOWER**
* **Optimization Modeling:** **YALMIP**
* **Solver:** **Gurobi Optimizer 10.0.2**

---

## 📊 Base System Data and Core Assumptions

### 1. Base Dataset

The basic system data is obtained from the **MATPOWER package**, utilizing the **IEEE 30-bus system** dataset.

### 2. Core Assumptions

The main modifications and core assumptions applied to the original IEEE 30-bus system are as follows:

* **Inverter-Based Resource Addition:** One IBR unit was added to the system to model a virtual inertia provider.
* **Demand and Cost Functions:** Load demand and the existing **generator cost functions** were kept **identical** to the MATPOWER base data.

### 3. Detailed Assumptions and Parameters

For **more detailed assumptions and data**, please refer to the separate file: **`ASSUMPTIONS.md`**.

