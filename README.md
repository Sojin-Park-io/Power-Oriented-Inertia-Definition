# Power-Oriented-Inertia-Definition

Supplementary material containing the complete mathematical formulation and simulation parameters for the paper:

**"Definition and Formulation of Inertia Service Incorporating Inverter-Based Resources"**

> **Status:** Submitted to the *IEEE Transactions on Power Systems* for peer review.

This repository serves as an archive for the mathematical models and simulation assumptions used to define and formalize inertia services incorporating inverter-based resources (IBR), as presented in the research.

---

## 📂 Repository Contents

This repository provides the theoretical backbone and data parameters of the research:

* **📄 Complete Formulation (PDF)**
    * File: **[`Formulation.pdf`](./Formulation.pdf)**
    * Description: Contains the detailed mathematical derivations and optimization problem formulations proposed in the paper.

* **📝 SCED Simulation Data & Assumptions**
    * Files: **[`ASSUMPTIONS_IEEE30.md`](./SCED%20Simulation/ASSUMPTIONS_IEEE30.md)** for the IEEE 30-bus system, **[`ASSUMPTIONS_IEEE118.md`](./SCED%20Simulation/ASSUMPTIONS_IEEE118.md)** for the IEEE 118-bus system
    * Description: Lists all specific simulation parameters, including grid modifications, generator characteristics, and cost function settings used to generate the results.

* **📝 TDS Simulation Data & Assumptions**
    * File: **[`ASSUMPTIONS.md`](./TDS%20Simulation/ASSUMPTIONS.md)** 
    * Description: Lists all specific simulation parameters, including governor time constants and recovery behavior settings.

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

## 🤝 Acknowledgements

This research and repository would not have been possible without the invaluable guidance and contribution of my co-authors.

I would like to express my deepest gratitude to **Prof. Ross Baldick** and **Prof. Hunyoung Shin** for their expert insight, supervision, and dedicated work on the paper, *"Definition and Formulation of Inertia Service Incorporating Inverter-Based Resources."* Their collaboration was essential to the successful completion of this project.
