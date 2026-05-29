## Detailed Assumptions for the Time-Domain Simulation (TDS)

This document provides a detailed breakdown of the simulation setup, model parameters, and scenario configurations used in the time-domain frequency simulations. The TDS uses dispatch results from the IEEE 118-bus SCED as its input.

---

### 1. Input from SCED

The TDS is driven by the pre-computed SCED dispatch results. The following quantities are loaded directly:

| Description |
| :--- |
| Per-Synchronous Generator (SG) Scheduled Output: $P_i^{\rm E}$, $P_i^{\rm In}$, $P_i^{\rm PFRr}$, $P_i^{\rm PFRd}$ |
| Per-Inverter-Based Resource (IBR) Scheduled Output: $P_j^{\rm E}$, $P_j^{\rm In}$, $P_j^{\rm PFRr}$, $P_j^{\rm PFRd}$ |
| Largest Contingency Size: $\Delta P^{\rm L}$ [MW] |
| IBR physical limits (capacity, SoC bounds): $P_j^{\min}$, $P_j^{\max}$, $E_j^{\min}$, $E_j^{\max}$, $E_{j0}$, $\eta$ |
| Frequency security limits used in SCED: `RoCoF_lim`, `Nadir_lim`, `QSS_lim` |
| IBR type indices (Type A / Type B): `iVI_only`, `iPFR_all` |

All TDS runs use Case 1 ($\alpha = 1.0$, $\beta = 0.0$) as the baseline SCED dispatch.

---

### 2. System Dynamics Model

All SGs are assumed to swing coherently, reducing the multi-machine system to a single equivalent swing equation:

$M_{\rm eff} \cdot \dot{f}(t) = -\Delta P^L + \sum_k P_k^E(t) + \sum_k P_k^{{\rm PFR}}(t) + \sum_{j,capped} P_j^{\rm In}(t)$

where $M_{\rm eff}$ is the effective inertia (SG + active IBR virtual inertia), computed iteratively to respect per-IBR headroom and footroom constraints.

#### 2a. Governor Model (SG)

SG primary frequency response is modeled as a first-order lag:

$$\dot{P}_{i}^{\rm PFR}(t) = \frac{P_{\rm sp,i}(t) - P_{i}^{\rm PFR}(t)}{T_{\rm SG}}$$

| Parameter | Value | Unit |
| :--- | :--- | :--- |
| **Governor time constant ($T_{\rm SG}$)** | $2.0$ | s |
| **Frequency deadband ($f_{\rm db}$)** | $0.036$ | Hz |

#### 2b. IBR Model

IBR PFR response is modeled as a first-order lag:

$$\dot{P}_{j}^{\rm PFR}(t) = \frac{P_{\rm sp,j}(t) - P_{j}^{\rm PFR}(t)}{T_{\rm IBR}}$$

| Parameter | Value | Unit |
| :--- | :--- | :--- |
| **IBR response time constant ($T_{\rm IBR}$)** | $0.1$ | s |
| **Round-trip efficiency ($\eta$)** | $0.9$ | — |

The setpoint $P_{\rm sp,j}$ is computed from the droop characteristic, clipped at the scheduled PFR capacity $P_j^{\rm PFR,r}$.

#### 2c. Virtual Inertia Response

The IBR virtual inertia output $P_{\rm VI,k}$ is determined in each time step, respecting:

- **Headroom:** $P_j^{\rm In} \leq \min(P_j^{\rm In}\; P_j^{\max} - P_j^{\rm E} - P_j^{\rm PFR}(t))$
- **Footroom:** $P_j^{\rm In} \geq P_j^{\rm In,low}$ (strategy-dependent)
- **SoC check:** $P_j^{\rm In} = 0$ if $E_j(t) \leq \underline{E}_j$

---

### 3. Simulation Setup

| Parameter | Value | Unit |
| :--- | :--- | :--- |
| **Simulation duration ($T_{\rm sim}$)** | $35$ | s |
| **ODE solver** | `ode45` | — |
| **Relative tolerance** | $10^{-7}$ | — |
| **Absolute tolerance** | $10^{-9}$ | — |
| **Maximum ODE step size** | $0.02$ | s |
| **Output time step** | $0.005$ | s |

---

### 4. Simulation Cases

Two TDS simulations are conducted, each examining a different reliability aspect of the P-oriented definition.

#### 4a. VI Recovery Strategies — Successive Contingency

This simulation compares three VI recovery strategies following a successive event:

| Event | Time | Description |
| :--- | :--- | :--- |
| **Event 1** | $t = 0 \text{s}$ | Largest generator trip ($\Delta P^L = 600 \text{MW}$) |
| **Event 2** | $t = 15 \text{s}$ | Second generator trip ($\Delta P^{L2} = 300 \text{MW}$) |

The three strategies compared are:

| Strategy | Description |
| :--- | :--- |
| **RoCoF Recovery** | $P_j^{\rm In} = -(M_{\rm In}/\eta)\cdot {\Delta f'}$ |
| **Constant Recovery (QSS)** | After each nadir, IBR absorbs a constant $P_j^{\rm const}$ for a fixed window ($t_{\rm delay} = 4 \text{s}$, $t_{\rm rec} = 4 \text{s}$) to recharge, then $P_j^{\rm In} = 0$ |
| **Discharge-only VI** | $P_j^{\rm In} \geq 0$ enforced (no absorption); IBR responds only when frequency is falling |

#### 4b. Symmetrical vs. Non-Symmetrical Footroom — Successive Event

This simulation compares two footroom allocation strategies under a successive gen-load trip:

| Event | Time | Description |
| :--- | :--- | :--- |
| **Event 1** | $t = 0 \text{s}$ | Generator trip ($\Delta P^{L1} = +600 \text{MW}$) |
| **Event 2** | $t = 15 \text{s}$ | Load trip ($\Delta P^{L2} = -500 \text{MW}$) |

The two footroom strategies compared are:

| Strategy | Footroom | Description |
| :--- | :--- | :--- |
| **Symmetrical** | $P_j^{\rm In} / \eta$ | Full absorption headroom; IBR can recover all discharged energy |
| **Non-Symmetrical** | $P_j^{\rm In} / (4\eta)$ | Reduced footroom (25% of symmetrical); clips VI absorption on frequency recovery |

When footroom is hit, the capped IBR's $M_j$ is excluded from the effective inertia, reducing effective absorption inertia and amplifying frequency overshoot.


## Time-Domain Simulation (TDS) Algorithm

### Overview

The TDS integrates four coupled dynamical systems simultaneously via MATLAB's `ode45` (adaptive step-size Runge-Kutta):

| State Variable | Description |
| :--- | :--- |
| $\Delta f(t)$ | COI frequency deviation |
| $P_i^{\text{Gov}}(t)$ | SG governor output (per SG) |
| $P_j^{\text{PFR}}(t)$ | IBR PFR output (per IBR) |
| $E_j^{\text{soc}}(t)$ | IBR SoC (per IBR) |

State vector dimension: $n_S = 1 + n_{SG} + 2n_{IBR}$

---

### Algorithm

```
INPUT:  SCED dispatch {P_k^E, P_k^In, P_k^PFRr, P_k^PFRd}
        Contingency size ΔP^L
        IBR parameters {P_j^min, P_j^max, η_j, E_j0, E_j^min, E_j^max}
        T_sg, T_ibr, f_db, T_sim

OUTPUT: Δf(t), {P_j^In(t)}, {P_j^PFR(t)}

─────────────────────────────────────────────────────────
INITIALIZATION
─────────────────────────────────────────────────────────
1. Set initial state:
   y₀ ← [Δf=0, {P_i^Gov=0}, {P_j^PFR=0}, {E_j^soc=0.5·E_j0}]

─────────────────────────────────────────────────────────
ODE INTEGRATION  (ode45 calls F(t,y) at each step)
─────────────────────────────────────────────────────────
2. Integrate ẏ = F(t, y) over [0, T_sim]

   ── F(t, y): Right-Hand Side ──────────────────────────

   3. Extract states:
      Δf, {P_i^Gov}, {P_j^PFR}, {E_j^soc} ← y

   4. SoC eligibility:
      s_j ← 1{ E_j^soc > E_j^min }          ▷ 0 if depleted

   5. Compute VI bounds:
      P̄_j^VI,+  ← max( min(P_j^In, P̄_j - P_j^ess - P_j^PFR), 0 )
      P_j^VI,-   ← min( -(P_j^ess + P_j^PFR - P_j^min), 0 )

   ── SWING EQUATION: Resolve Δf' and P_j^VI ──────────
   (Fixed-point iteration with saturation sweep)

   6. Initialize: capped_j ← false,  P_j^VI,cap ← 0

   7. for iter = 1 to n_IBR + 2:

      a. free_j ← s_j AND NOT capped_j

      b. rhs ← -ΔP^L + Σ P_i^Gov + Σ(P_j^PFR · s_j) + Σ P_j^VI,cap

      c. if rhs ≤ 0:                          ▷ FDP (frequency decline)
            M_eff ← M_sg + Σ(M_j^vi · free_j)
            Δf'   ← rhs / M_eff
            P_j^VI,trial ← M_j^vi · free_j · (-Δf')

         else:                                ▷ FRP (frequency recovery)
            M_eff ← M_sg + Σ(M_j^vi/η_j · free_j)
            Δf'   ← rhs / M_eff
            P_j^VI,trial ← -(M_j^vi/η_j) · free_j · Δf'

      d. for each j:
            if P_j^VI,trial > P̄_j^VI,+:      ▷ Headroom clipping
               capped_j ← true
               P_j^VI,cap ← P̄_j^VI,+

            elif P_j^VI,trial < P_j^VI,-:     ▷ Footroom clipping
               capped_j ← true
               P_j^VI,cap ← P_j^VI,-

            else: break                       ▷ No violations → converged

   8. Assign final VI output:
      P_j^VI ← P_j^VI,trial   if free_j
             ← P_j^VI,cap     if capped_j
             ← 0              if NOT s_j

   ── PFR SETPOINTS (first-order lag with deadband) ────

   9.  Δf_exc ← max(-Δf - f_db, 0)

   10. P_i^Gov,sp ← min(K_i^d · Δf_exc,  P̄_i^PFRr)

   11. P_j^PFR,sp ← s_j · min(K_j^d · Δf_exc,
                               P̄_j^PFRr,
                               max(P̄_j - P_j^ess - P_j^VI, 0))

   ── ASSEMBLE ẏ ───────────────────────────────────────

   12. dΔf/dt      ← Δf'
       dP_i^Gov/dt ← (P_i^Gov,sp - P_i^Gov)  / T_sg
       dP_j^PFR/dt ← (P_j^PFR,sp - P_j^PFR) / T_ibr
       dE_j^soc/dt ← -(P_j^ess + P_j^VI + s_j·P_j^PFR) / 3600

   13. return ẏ
```

---

### Key Design Notes

| Feature | Description |
| :--- | :--- |
| **VI resolution** | Algebraic fixed-point at each integration step; no separate ODE for VI |
| **FDP/FRP branching** | Different effective inertia depending on sign of `rhs`: discharge (`M_vi`) vs. absorption (`M_vi/η`) |
| **Saturation sweep** | Monotone iteration — once an IBR is capped, it stays capped; converges within `n_IBR + 2` iterations |
| **Headroom clipping** | `P_j^VI` capped at `P̄_j - P_j^ess - P_j^PFR` |
| **Footroom clipping** | `P_j^VI` floored at `-(P_j^ess + P_j^PFR - P_j^min)` |
| **SoC gate** | IBR excluded from both VI and PFR if `E_j^soc ≤ E_j^min` |
