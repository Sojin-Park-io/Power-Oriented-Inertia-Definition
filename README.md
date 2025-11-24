# Power-Oriented-Inertia-Definition

Supplementary material containing the complete mathematical formulation and simulation parameters for the paper:

**"Definition and Formulation of Inertia Service Incorporating Inverter-Based Resources"**

> **Status:** Submitted to the *IEEE Transactions on Power Systems* for peer review.

This repository serves as an archive for the mathematical models and simulation assumptions used to define and formalize inertia services incorporating inverter-based resources (IBR), as presented in the research.

---

## 📂 Repository Contents

This repository provides the theoretical backbone and data parameters of the research:

* **📄 Complete Formulation (PDF)**
    * File: **[`Formulation.pdf`}(./Formulation.pdf)**
    * Description: Contains the detailed mathematical derivations and optimization problem formulations proposed in the paper.

* **📝 Simulation Data & Assumptions**
    * File: **[`ASSUMPTIONS.md`](./ASSUMPTIONS.md)**
    * Description: Lists all specific simulation parameters, including grid modifications, generator constants, and cost function settings used to generate the results.

---

## ⚙️ Simulation Environment

The simulation results and optimization solutions discussed in the paper (and the provided formulation) were generated using the following environment. This information is provided to clarify the computational context of the research.

### 💻 Hardware
* **CPU:** Apple M3 chip (4.05 GHz clock rate)
* **RAM:** 16 GB Unified Memory

### 🛠️ Software & Solvers
* **Main Environment:** MATLAB (R2025a)
* **Modeling Framework:** YALMIP
* **Solver:** Gurobi Optimizer 10.0.2
* **Power System Data Source:** MATPOWER

---

## 📊 Base System Overview

### 1. Base Dataset
The study utilizes the **IEEE 30-bus system** dataset obtained from the **MATPOWER** package as the baseline.

### 2. Key Modifications
To analyze the impact of inertia services, the following modifications were applied to the standard IEEE 30-bus system:

* **Synchronous Inertia (SI) Scenarios:** The inertia constant ($H$) of synchronous generators was explicitly adjusted to model **Low-SI** and **High-SI** grid scenarios.
* **IBR Integration:** A virtual inertia provider (IBR unit) was integrated into the system model.
* **Parameters:** Load demand and generator cost functions remain **identical** to the original MATPOWER specifications.

> **Note:** For the specific numerical values of these modifications (e.g., exact $H$ values, IBR capacity), please refer to the **[`ASSUMPTIONS.md`](./ASSUMPTIONS.md)** file.

---

## 🤝 Acknowledgements

This research and repository would not have been possible without the invaluable guidance and contribution of my co-authors.

I would like to express my deepest gratitude to **Prof. Ross Baldick** and **Prof. Hunyoung Shin** for their expert insight, supervision, and dedicated work on the paper, *"Definition and Formulation of Inertia Service Incorporating Inverter-Based Resources."* Their collaboration was essential to the successful completion of this project.
