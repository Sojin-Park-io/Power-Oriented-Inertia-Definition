## Detailed Assumptions and Parameter Adjustments

This document provides a detailed breakdown of the adjustments made to the base IEEE 30-bus system data, particularly focusing on the parameters related to inertia, the Inverter-Based Resource (IBR), and operational constraints.

---

### 1. Synchronous Inertia (SI) Level Adjustment

The base system's synchronous inertia was explicitly modified to define the two primary grid scenarios discussed in the paper. This adjustment was achieved by setting the **Inertia Constant ($H$)** for the **six synchronous generators (SGs)**.

The two distinct SI levels are defined as follows:

| Scenario | Synchronous Inertia Constant ($H$) |
| :--- | :--- |
| **Low-SI Grid** | $H = 1.0 \, \text{s}$ for all six SGs |
| **High-SI Grid** | $H = 5.0 \, \text{s}$ for all six SGs |

---

### 2. Inverter-Based Resource (IBR) Integration (Virtual Inertia Provider)

One IBR unit, which **incorporates an Energy Storage System (ESS)**, is added to the system to model a virtual inertia provider. The IBR is assumed to be capable of providing **both positive and negative inertia responses** (like those achieved through VSG/VSM control). The detailed characteristics of the integrated ESS are as follows:

| Parameter | Value | Unit | Reference |
| :--- | :--- | :--- | :--- |
| **ESS Capacity** | $100 / 100$ | MW / MWh | - |
| **Round-Trip Efficiency ($\eta$)** | $0.9$ | - | - |
| **Min State of Charge ($\text{SoC}_{\min}$)** | $0.2$ | p.u. | - |
| **Max State of Charge ($\text{SoC}_{\max}$)** | $0.8$ | p.u. | - |
| **Operational Restriction** | Discharge-only operation | - | - |
| **Variable Energy Cost** | $3$ | \$/MWh | [1] |
| **Inertia Constant ($H$) Range** | $0 - 50$ | s | [2] |

The variable energy cost of \$3/MWh is applied to the IBR for both energy dispatch and the loss term, considering the degradation cost of the ESS.

---

### 3. System Load and Operational Parameters

The system operation is based on the Economic Dispatch (ED) problem following a prior UC solution where all generators are online.

* **System Demand:** $189.2 \, \text{MW}$ (approx. 56% of total capacity).
* **Droop Constant ($R$):** $5\%$ for all resources.
* **Dispatch Interval ($\Delta T$):** $5$ minutes.
* **Required Time for Fully-Activated PFR ($T^{\rm PFR}$):** $6$ seconds [3].
* **Nominal Frequency ($f_0$):** $60 \, \text{Hz}$.

### 4. Frequency Constraint Limits

The optimization formulation uses the following maximum allowable limits for frequency deviations:

| Constraint | Limit Value | Unit |
| :--- | :--- | :--- |
| **Maximum Allowable RoCoF ($\overline{\Delta f'}$)** | $1.0$ | Hz/s |
| **Maximum Allowable Frequency Deviation at Nadir ($\overline{\Delta f}_{\rm Nad}$)** | $0.8$ | Hz |
| **Maximum Allowable Frequency Deviation at QSS ($\overline{\Delta f}_{\rm QSS}$)** | $0.5$ | Hz |

---

## 📚 References

* **[1]** P. Kushwaha, V. Prakash, and R. Bhakar, "A novel framework to assess synthetic inertia \& primary frequency response support from energy storagy systems," \textit{Sustainable Energy, Grids and Networks}, vol. 34, 6 2023.
* **[2]** AEMO, "Hornsdale power reserve virtual machine mode testing summary report,” AEMO, Tech. Rep., Aug. 2023, issued August 2023.
* **[3]** P.S AEMO, \textit{Market Ancillary Services Specification (MASS), Version 8.2}, AEMO, June 2024, effective 3 June 2024.
