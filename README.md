# CSTR Performance Simulation — First-Order Reaction Kinetics

> **Independent Chemical Engineering Project**
> Reaction Engineering · Python Numerical Modeling · Process Design
> *Chemistry and Engineering of Organic Compounds | Petrochemistry and Carbochemistry*

---

## What This Project Does

Simulates a **Continuous Stirred Tank Reactor (CSTR)** for a first-order irreversible liquid-phase reaction **A → B**, analyzes how conversion responds to residence time and rate constant, validates the analytical solution through numerical iteration, and produces engineering-grade visualizations.

---

## Engineering Background

For a CSTR at steady state, the material balance gives the design equation:

```
X = (k · τ) / (1 + k · τ)
```

| Symbol | Meaning | Unit |
|--------|---------|------|
| X | Fractional conversion of A | — |
| k | First-order rate constant | min⁻¹ |
| τ | Residence time (V/F) | min |
| Da | Damköhler number = k·τ | — |

The **Damköhler number** is dimensionless — it universally characterizes CSTR performance regardless of individual k and τ. Da = 1 gives 50% conversion. Da = 9 gives 90%.

---

## Results

### Conversion Table (k = 0.5 min⁻¹)

| τ (min) | Conversion X | Remaining C_A/C_A0 |
|---------|-------------|---------------------|
| 1.0 | 0.3333 | 0.6667 |
| 2.0 | 0.5000 | 0.5000 |
| 5.0 | 0.7143 | 0.2857 |
| 10.0 | 0.8333 | 0.1667 |
| 20.0 | 0.9091 | 0.0909 |

### Numerical Validation

Analytical and iterative solutions converge within **< 1×10⁻⁸** error — confirming correctness of both methods.

### Design Example

- k = 0.5 min⁻¹, τ = 5 min → **X = 71.4%**
- To achieve **80% conversion**: required volume = **160 L** (at F_A0 = 10 mol/min, C_A0 = 1 mol/L)

---

## Plots Generated

| Plot | Description |
|------|-------------|
| Conversion vs τ | Parametric curves for 5 different k values |
| Universal CSTR curve | Conversion vs Damköhler number |
| Volume vs conversion | Reactor sizing curve |

---

## How to Run

```bash
pip install numpy matplotlib
python cstr_simulation.py
```

Output: terminal report + `cstr_results.png`

---

## Relevance to Petrochemistry

CSTRs are widely used in **organic compound synthesis and polymer production** — for example, continuous polymerization reactors, alkylation units, and esterification processes. Understanding conversion-residence time relationships is fundamental to sizing these reactors in petrochemical plant design.

---

## Skills Demonstrated

- CSTR design equation derivation and application
- Damköhler number analysis
- Numerical iteration with convergence validation
- Python scientific computing (NumPy, Matplotlib)
- Reactor volume sizing for target conversion

---

*Zidanur Rahman | Independent Chemical Engineering Study | 2025*
*Prepared as part of Romanian Government Scholarship application*
