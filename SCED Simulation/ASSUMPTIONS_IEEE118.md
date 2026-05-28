## Detailed Assumptions for the IEEE 118-Bus System and Parameter Adjustments

This document provides a detailed breakdown of the adjustments made to the IEEE 118-bus system data, particularly focusing on the parameters related to inertia, the Inverter-Based Resources (IBRs), and operational constraints.

---

### 1. Synchronous Generator (SG) Configuration

The IEEE 118-bus base system includes 54 generators total. The last 12 are re-designated as IBRs, resulting in 42 synchronous generators (SGs). The generator with the maximum capacity is set to $600 \, \text{MW}$ and designated as the reference contingency unit; all remaining SGs are scaled down to adjust the IBR penetration share to the target level (35%). The modified capacity values for each unit are provided in the accompanying file [GenData_IEEE118.xlsx](./GenData_IEEE118.xlsx).


The inertia constant is uniformly assigned to all SGs as follows:

| Parameter | Value | Unit |
| :--- | :--- | :--- |
| **Inertia Constant ($H$)** | $5.0$ | s |

---

### 2. IBR Integration

Twelve IBR units, comprising two distinct types, are integrated into the system to serve as virtual inertia providers and ancillary service contributors.

#### 2a. IBR Type Classification

| Type | Units | Service Capability | $H_{\text{IBR,max}}$ | PFR Capability |
| :--- | :--- | :--- | :--- | :--- |
| **Type A** | IBR 1–3 (3 units) | Energy + Inertia only | $50 \text{s}$ | None |
| **Type B** | IBR 4–12 (9 units) | Energy + Inertia + PFR | $10 \text{s}$ | Available |

* **Type A** units carry a large VI headroom ($H_{\text{max}} = 50 \text{s}$) in exchange for no PFR obligation, acting as dedicated inertia providers.
* **Type B** units provide all ancillary services with a limited VI capability ($\(H_{\text{max}} = 10 \text{s}$)

#### 2b. Common IBR Parameters

The following parameters apply to **all 12 IBR units**:

| Parameter | Value | Unit | Reference |
| :--- | :--- | :--- | :--- |
| **Energy Capacity ($E_{\text{cap}}$)** | $P_{\max} \times 1$ | MW / MWh (1-hour rated) | — |
| **Round-Trip Efficiency ($\eta$)** | $0.9$ | — | — |
| **Min State of Charge ($\text{SoC}_{\min}$)** | $0.2$ | p.u. | — |
| **Max State of Charge ($\text{SoC}_{\max}$)** | $0.8$ | p.u. | — |
| **Initial State of Charge ($\text{SoC}_{0}$)** | $0.5$ | p.u. | — |
| **Minimum Inertia Constant ($H_{\text{min}}$)** | $0$ | s | — |
| **Variable Energy Cost** | $3$ | \$/MWh | [1] |

The variable energy cost of \$3/MWh is applied to both the energy dispatch and the efficiency-loss term to account for ESS degradation costs.

---

### 3. System Capacity Configuration

| Parameter | Value |
| :--- | :--- |
| **Total System Load** | $4{,}831.0 \text{MW}$ |
| **Net Load** | $4{,}242.0 \text{MW}$ |
| **RE Output (non-dispatchable, CF = 100%)** | $589.0 \text{MW}$ |
| **IBR dispatchable share** | 30% |
| **IBR total penetration** | 35% |
| **Fixed contingency size** | $600 \text{MW}$ |

The net load served by the SCED equals the total IEEE 118-bus system demand (approximately $4,242 \text{MW}$), after subtracting the non-dispatchable RE output, which is treated as a negative load.

---

### 4. System Load and Operational Parameters

The system operation is based on the Economic Dispatch (ED) problem following a prior Unit Commitment (UC) solution where all generators are online.

* **Droop Constant ($R$):** $5\%$.
* **Dispatch Interval ($\Delta T$):** $5$ minutes.
* **Required Time for Fully-Activated PFR ($T^{\rm PFR}$):** $6$ seconds [3].
* **Nominal Frequency ($f_0$):** $60 \text{Hz}$.

---

### 5. Frequency Constraint Limits

The optimization formulation uses the following frequency security limits.

| Constraint | Limit Value | Unit |
| :--- | :--- | :--- |
| **Maximum Allowable RoCoF ($\overline{\Delta f'}$)** | $0.5$ | Hz/s |
| **Maximum Allowable Frequency Deviation at Nadir ($\overline{\Delta f}_{\rm Nad}$)** | $0.5$ | Hz |
| **Maximum Allowable Frequency Deviation at QSS ($\overline{\Delta f}_{\rm QSS}$)** | $0.3$ | Hz |

---

### 6. Simulation Cases

#### 6a. Largest Contingency Formulation Cases

Two approaches to contingency modeling are compared, with $\alpha = 1.0$ and $\beta = 0.0$ fixed:

| Case | Description |
| :--- | :--- |
| **Variable Contingency** | The largest contingency $\Delta P^{L}$ is a decision variable; SCED endogenously determines the binding contingency unit |
| **Fixed Contingency** | $\Delta P^{L} = \max(P_{\max}) = 600 \text{MW}$ — imposes the worst-case contingency regardless of dispatch |

#### 6b. $\alpha$/$\beta$ Sensitivity Analysis

Under the fixed contingency setting ($\Delta P^{L} = 600 \text{MW}$), the following two cases are compared to examine the effect of $\alpha$ and $\beta$ on headroom and footroom constraints:

| Case | $\alpha$ | $\beta$ | Description |
| :--- | :--- | :--- | :--- |
| **Case 1** | $1.0$ | $0.0$ | Conservative — full headroom reserved for both inertia and PFR services; ramp PFR contributes nothing to footroom relief |
| **Case 2** | $0.5$ | $0.5$ | Flexible — partial headroom reserved for both inertia and PFR services; ramp PFR partially relieves the footroom requirement |


---

## 📚 References

* **[1]** P. Kushwaha, V. Prakash, and R. Bhakar, "A novel framework to assess synthetic inertia \& primary frequency response support from energy storage systems," *Sustainable Energy, Grids and Networks*, vol. 34, Jun. 2023.
* **[2]** AEMO, "Hornsdale power reserve virtual machine mode testing summary report," AEMO, Tech. Rep., Aug. 2023, issued August 2023.
* **[3]** AEMO, *Market Ancillary Services Specification (MASS), Version 8.2*, AEMO, June 2024, effective 3 June 2024.
