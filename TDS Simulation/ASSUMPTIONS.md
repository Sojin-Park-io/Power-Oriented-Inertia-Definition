# Time-Domain Simulation (TDS) — Detailed Documentation

This document provides a detailed breakdown of the simulation setup, model parameters, scenario configurations, and algorithm used in the time-domain frequency simulations. The TDS uses dispatch results from the IEEE 118-bus SCED as its input.

---

## 1. Input from SCED

The TDS is driven by the pre-computed SCED dispatch results. The following quantities are loaded directly:

| Description |
| :--- |
| Per-SG scheduled output: $P_i^{\rm E}$, $P_i^{\rm In}$, $P_i^{\rm PFR,r}$, $P_i^{\rm PFR,d}$ |
| Per-IBR scheduled output: $P_j^{\rm E}$, $P_j^{\rm In}$, $P_j^{\rm PFR,r}$, $P_j^{\rm PFR,d}$ |
| Largest contingency size: $\Delta P^{\rm L}$ |
| IBR physical limits: $P_j^{\min}$, $P_j^{\max}$, $E_j^{\min}$, $E_j^{\max}$, $E_{j0}$, $\eta_j$ |
| Frequency security limits used in SCED: `RoCoF_lim`, `Nadir_lim`, `QSS_lim` |
| IBR type indices (Type A / Type B): `iVI_only`, `iPFR_all` |

All TDS runs use Case 1 ($\alpha = 1.0$, $\beta = 0.0$) as the baseline SCED dispatch.

---
| **SCED vs. TDS variables** | $P_j^{\rm In}$, $P_j^{\rm E}$, $P_j^{\rm PFR,r}$ are SCED scalars (no (t)); $P_j^{\rm In}(t)$, $P_j^{\rm PFR}(t)$, $\Delta f(t)$ are TDS time functions |

## 2. System Dynamics Model

All resources are assumed to swing coherently, reducing the multi-machine system to a single equivalent swing equation:

$$M_{\rm eff}(t) \cdot \Delta f'(t) = -\Delta P^{\rm L} + \sum_i P_i^{\rm PFR}(t) + \sum_j P_j^{\rm PFR}(t) \cdot s_j(t) + \sum_j P_j^{\rm In}(t)$$

where $M_{\rm eff}(t)$ is the effective inertia (SG + active IBR virtual inertia), computed iteratively at each time step to respect per-IBR headroom and footroom constraints.

### 2a. Governor Model (SG)

SG primary frequency response is modeled as a first-order lag:

$$\dot{P}_{i}^{\rm PFR}(t) = \frac{P_{i}^{\rm PFR,sp}(t) - P_{i}^{\rm PFR}(t)}{T_{\rm SG}}$$

| Parameter | Value | Unit |
| :--- | :--- | :--- |
| Governor time constant ($T_{\rm SG}$) | $2.0$ | s |
| Frequency deadband ($f_{\rm db}$) | $0.036$ | Hz |

The setpoint $P_i^{\rm PFR,sp}(t)$ is computed from the droop characteristic, clipped at the scheduled PFR capacity $P_i^{\rm PFR,r}$:

$$P_i^{\rm PFR,sp}(t) = \min \left(K_i^d \cdot \Delta f_{\rm exc}(t),\ P_i^{\rm PFR,r}\right)$$

$$\Delta f_{\rm exc}(t) = \max(-\Delta f(t) - f_{\rm db},\ 0)$$

### 2b. IBR PFR Model

IBR PFR response is modeled as a first-order lag:

$$\dot{P}_{j}^{\rm PFR}(t) = \frac{P_{j}^{\rm PFR,sp}(t) - P_{j}^{\rm PFR}(t)}{T_{\rm IBR}}$$

| Parameter | Value | Unit |
| :--- | :--- | :--- |
| IBR response time constant ($T_{\rm IBR}$) | $0.1$ | s |
| Round-trip efficiency ($\eta_j$) | $0.9$ | — |

The setpoint $P_j^{\rm PFR,sp}(t)$ is clipped at the scheduled PFR capacity and the residual headroom after VI:

$$P_j^{\rm PFR,sp}(t) = \min \left(K_j^d \cdot \Delta f_{\rm exc}(t),\ P_j^{\rm PFR,r},\ \max(P_j^{\max} - P_j^{\rm E} - P_j^{\rm In}(t),\ 0)\right)$$

### 2c. Virtual Inertia Response

The IBR virtual inertia output $P_j^{\rm In}(t)$ is determined at each time step, respecting:

- **Headroom:** $P_j^{\rm In}(t) \leq \max\left(\min(P_j^{\rm In},\ P_j^{\max} - P_j^{\rm E} - P_j^{\rm PFR}(t)),\ 0\right)$
- **Footroom:** $P_j^{\rm In}(t) \geq \min\left(-(P_j^{\rm E} + P_j^{\rm PFR}(t) - P_j^{\min}),\ 0\right)$
- **SoC check:** $P_j^{\rm In}(t) = 0$ if $E_j(t) \leq E_j^{\min}$

### 2d. SoC Dynamics

$$\dot{E}_j(t) = -\frac{1}{3600}\left(P_j^{\rm E} + P_j^{\rm In}(t) + P_j^{\rm PFR}(t)\right)$$

---

## 3. Simulation Setup

| Parameter | Value | Unit |
| :--- | :--- | :--- |
| Simulation duration ($T_{\rm sim}$) | $35$ | s |
| ODE solver | `ode45` | — |
| Relative tolerance | $10^{-7}$ | — |
| Absolute tolerance | $10^{-9}$ | — |
| Maximum ODE step size | $0.02$ | s |
| Output time step | $0.005$ | s |

---

## 4. Simulation Cases

Two TDS simulations are conducted, each examining a different reliability aspect of the P-oriented definition.

### 4a. VI Recovery Strategies — Successive Contingency

This simulation compares three VI recovery strategies following successive generation trip events:

| Event | Time | Description |
| :--- | :--- | :--- |
| Event 1 | $t = 0$ s | Largest generator trip ($\Delta P^{\rm L} = 600$ MW) |
| Event 2 | $t = 15$ s | Second generator trip ($\Delta P^{\rm L2} = 300$ MW) |

| Strategy | Description |
| :--- | :--- |
| **RoCoF Recovery** | $P_j^{\rm In}(t) = -(M_j^{\rm In}/\eta_j) \cdot \Delta f'(t)$ in both FDP and FRP |
| **Constant Recovery** | After each nadir, IBR absorbs a constant $P_j^{\rm const}$ for a fixed window ($t_{\rm delay} = 4$ s, $t_{\rm rec} = 4$ s), then $P_j^{\rm In}(t) = 0$ |
| **No Recovery** | $P_j^{\rm In}(t) \geq 0$ enforced; IBR responds only when frequency is declining |

### 4b. Symmetrical vs. Non-Symmetrical Footroom — Successive Event

This simulation compares two footroom allocation strategies under a successive generation and load trip:

| Event | Time | Description |
| :--- | :--- | :--- |
| Event 1 | $t = 0$ s | Generator trip ($\Delta P^{\rm L1} = +600$ MW) |
| Event 2 | $t = 15$ s | Load trip ($\Delta P^{\rm L2} = -500$ MW) |

| Strategy | Footroom | Description |
| :--- | :--- | :--- |
| **Symmetrical** | $P_j^{\rm In} / \eta_j$ | Full absorption capacity; IBR can recover all discharged energy |
| **Non-Symmetrical** | $P_j^{\rm In} / (4\eta_j)$ | Reduced footroom (25% of symmetrical); clips $P_j^{\rm In}(t)$ during FRP |

When footroom is hit, the capped IBR's $M_j^{\rm In}$ is excluded from $M_{\rm eff}(t)$, reducing effective absorption inertia and amplifying frequency overshoot.

---

## 5. TDS Algorithm

### State Vector

| Variable | Description | Dimension |
| :--- | :--- | :--- |
| $\Delta f(t)$ | COI frequency deviation | $1$ |
| $P_i^{\rm PFR}(t)$ | SG governor (PFR) output | $n_{SG}$ |
| $P_j^{\rm PFR}(t)$ | IBR PFR output | $n_{IBR}$ |
| $E_j(t)$ | IBR state of charge | $n_{IBR}$ |

Total dimension: $n_S = 1 + n_{SG} + 2n_{IBR}$

---

```
INPUT
  SCED dispatch : {P_k^E, P_k^In, P_k^PFRr, P_k^PFRd}   ← scalars, no (t)
  Contingency   : ΔP^L
  IBR params    : {P_j^min, P_j^max, η_j, E_j0, E_j^min}
  Dynamics      : T_SG, T_IBR, f_db, T_sim

OUTPUT
  Δf(t),  {P_j^In(t)},  {P_j^PFR(t)}                    ← time functions

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 1. INITIALIZATION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  y₀ ← [Δf=0,  {P_i^PFR=0},  {P_j^PFR=0},  {E_j=E_j0}]

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
STEP 2. ODE INTEGRATION  (ode45, adaptive RK)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Integrate ẏ = F(t, y) over [0, T_sim]
  F(t, y) is evaluated at each integration step:

  ┌─────────────────────────────────────────────┐
  │  F(t, y) : Right-Hand Side                  │
  └─────────────────────────────────────────────┘

  2-1. Extract states from y:
       Δf(t),  {P_i^PFR(t)},  {P_j^PFR(t)},  {E_j(t)}

  2-2. SoC eligibility:
       s_j(t) ← 1{ E_j(t) > E_j^min }          ▷ 0 if depleted

  2-3. Compute VI bounds at current time step:
       P_j^{In,+}(t) ← max( min(P_j^In, P_j^max - P_j^E - P_j^PFR(t)), 0 )
       P_j^{In,-}(t) ← min( -(P_j^E + P_j^PFR(t) - P_j^min), 0 )

  ┌─────────────────────────────────────────────┐
  │  Swing Equation                             │
  │  Resolve Δf'(t) and P_j^In(t)               │
  └─────────────────────────────────────────────┘

  2-4. Initialize:
       capped_j ← false,  P_j^{In,cap}(t) ← 0

  2-5. for iter = 1 to n_IBR + 2:

       free_j ← s_j(t) AND NOT capped_j

       rhs ← -ΔP^L + Σ_i P_i^PFR(t) + Σ_j (P_j^PFR(t) · s_j(t)) + Σ_j P_j^{In,cap}(t)

       if rhs ≤ 0:                      ▷ FDP (frequency decline)
         M_eff ← M_sg + Σ_j(M_j^In · free_j)
         Δf'(t)         ← rhs / M_eff
         P_j^{In,trial}(t) ← M_j^In · free_j · (-Δf'(t))

       else:                            ▷ FRP (frequency recovery)
         M_eff ← M_sg + Σ_j(M_j^In/η_j · free_j)
         Δf'(t)         ← rhs / M_eff
         P_j^{In,trial}(t) ← -(M_j^In/η_j) · free_j · Δf'(t)

       for each j:
         if P_j^{In,trial}(t) > P_j^{In,+}(t):   ▷ Headroom clipping
           capped_j       ← true
           P_j^{In,cap}(t) ← P_j^{In,+}(t)

         elif P_j^{In,trial}(t) < P_j^{In,-}(t):  ▷ Footroom clipping
           capped_j       ← true
           P_j^{In,cap}(t) ← P_j^{In,-}(t)

         else: break

  2-6. Assign final VI output:
       P_j^In(t) ← P_j^{In,trial}(t)  if free_j      (unconstrained)
                ← P_j^{In,cap}(t)     if capped_j    (at headroom/footroom limit)
                ← 0                   if NOT s_j(t)  (SoC depleted)

  ┌─────────────────────────────────────────────┐
  │  PFR Setpoints                              │
  │  First-order lag with deadband              │
  └─────────────────────────────────────────────┘

  2-7. Δf_exc(t) ← max(-Δf(t) - f_db,  0)

  2-8. P_i^{PFR,sp}(t) ← min(K_i^d · Δf_exc(t),  P_i^{PFR,r})

  2-9. P_j^{PFR,sp}(t) ← s_j(t) · min(K_j^d · Δf_exc(t), P_j^{PFR,r}, max(P_j^max - P_j^E - P_j^In(t), 0))

  ┌─────────────────────────────────────────────┐
  │  Assemble ẏ                                 │
  └─────────────────────────────────────────────┘

  2-10. dΔf/dt        ← Δf'(t)
        dP_i^PFR/dt   ← (P_i^{PFR,sp}(t) - P_i^PFR(t))  / T_SG
        dP_j^PFR/dt   ← (P_j^{PFR,sp}(t) - P_j^PFR(t))  / T_IBR
        dE_j/dt       ← -(P_j^E + P_j^In(t) + s_j(t)·P_j^PFR(t)) / 3600

  return ẏ
```

---



