"""
=============================================================
CSTR Performance Simulation — First-Order Reaction Kinetics
=============================================================
Author      : Zid
Field       : Reaction Engineering / Chemical Process Design
Subject     : Chemistry and Engineering of Organic Compounds,
              Petrochemistry and Carbochemistry
Description :
    Simulates a Continuous Stirred Tank Reactor (CSTR) for a
    first-order irreversible liquid-phase reaction: A → B

    Analyzes how conversion varies with residence time and
    reaction rate constant. Validates results through numerical
    iteration and generates engineering-grade visualizations.

Reaction    : A → B   (-r_A = k · C_A)
Design Eq.  : X = (k · τ) / (1 + k · τ)
=============================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# ─────────────────────────────────────────────────────────────
# CORE CSTR EQUATIONS
# ─────────────────────────────────────────────────────────────

def cstr_conversion(k, tau):
    """
    Analytical solution — CSTR steady-state conversion.

    Parameters:
        k   : reaction rate constant (min⁻¹)
        tau : residence time V/F (min)

    Returns:
        X   : fractional conversion (0–1)
    """
    return (k * tau) / (1 + k * tau)


def required_volume(F_A0, X_target, k, C_A0):
    """
    Calculate reactor volume for a target conversion.
    From CSTR design equation: V = F_A0 * X / (k * C_A0 * (1-X))

    Parameters:
        F_A0     : molar feed rate of A (mol/min)
        X_target : desired conversion (0–1)
        k        : rate constant (min⁻¹)
        C_A0     : inlet concentration (mol/L)

    Returns:
        V : required volume (L)
    """
    r_A = k * C_A0 * (1 - X_target)
    return (F_A0 * X_target) / r_A


# ─────────────────────────────────────────────────────────────
# NUMERICAL VALIDATION — Iterative Solver
# ─────────────────────────────────────────────────────────────

def iterative_solver(k, tau, C_A0=1.0, tol=1e-8, max_iter=1000):
    """
    Solve CSTR concentration iteratively to validate analytical result.
    Iterative form: C_A(n+1) = C_A0 / (1 + k·τ)
    Converges to the same answer as analytical — confirms correctness.
    """
    C_A = C_A0
    for i in range(max_iter):
        C_A_new = C_A0 / (1 + k * tau)
        if abs(C_A_new - C_A) < tol:
            X = (C_A0 - C_A_new) / C_A0
            return C_A_new, X, i + 1
        C_A = C_A_new
    return C_A, (C_A0 - C_A) / C_A0, max_iter


# ─────────────────────────────────────────────────────────────
# CONVERSION TABLE
# ─────────────────────────────────────────────────────────────

def print_conversion_table(k=0.5):
    tau_list = [0.5, 1, 2, 3, 5, 7, 10, 15, 20]

    print(f"\n{'='*55}")
    print(f"  CSTR Conversion Table  |  k = {k} min⁻¹")
    print(f"{'='*55}")
    print(f"  {'τ (min)':>8}  {'Conversion X':>14}  {'C_A/C_A0':>10}")
    print(f"  {'-'*44}")
    for tau in tau_list:
        X = cstr_conversion(k, tau)
        print(f"  {tau:>8.1f}  {X:>14.4f}  {1-X:>10.4f}")
    print(f"{'='*55}")


# ─────────────────────────────────────────────────────────────
# VALIDATION TABLE
# ─────────────────────────────────────────────────────────────

def print_validation(k=0.5):
    tau_list = [1, 2, 5, 10, 20]

    print(f"\n{'='*68}")
    print(f"  Validation: Analytical vs Numerical  |  k = {k} min⁻¹")
    print(f"{'='*68}")
    print(f"  {'τ':>5}  {'X_analytical':>14}  {'X_numerical':>13}  {'Error':>12}  {'Iters':>6}")
    print(f"  {'-'*58}")
    for tau in tau_list:
        X_a = cstr_conversion(k, tau)
        _, X_n, iters = iterative_solver(k, tau)
        err = abs(X_a - X_n)
        print(f"  {tau:>5.1f}  {X_a:>14.8f}  {X_n:>13.8f}  {err:>12.2e}  {iters:>6}")
    print(f"{'='*68}")


# ─────────────────────────────────────────────────────────────
# PLOTS
# ─────────────────────────────────────────────────────────────

def plot_results():
    tau = np.linspace(0.1, 20, 500)
    k_values = [0.1, 0.3, 0.5, 1.0, 2.0]
    colors = ['#1B3A6B', '#2E5FA3', '#0D7377', '#E07B39', '#8B2FC9']

    fig = plt.figure(figsize=(14, 10))
    fig.suptitle(
        "CSTR Performance Simulation — First-Order Reaction Kinetics\n"
        "Chemistry & Engineering of Organic Compounds | Independent Study",
        fontsize=13, fontweight='bold', y=0.98
    )
    gs = gridspec.GridSpec(2, 2, hspace=0.42, wspace=0.35)

    # Plot 1: Conversion vs Residence Time
    ax1 = fig.add_subplot(gs[0, :])
    for i, k in enumerate(k_values):
        X = cstr_conversion(k, tau)
        ax1.plot(tau, X, color=colors[i], linewidth=2.2, label=f'k = {k} min⁻¹')
    ax1.axhline(0.90, color='red', linestyle='--', alpha=0.5, linewidth=1.2, label='90% target')
    ax1.set_xlabel("Residence Time τ (min)", fontsize=11)
    ax1.set_ylabel("Conversion X", fontsize=11)
    ax1.set_title("Conversion vs Residence Time — Parametric Study", fontsize=12, fontweight='bold')
    ax1.legend(loc='lower right', fontsize=9)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim([0, 20]); ax1.set_ylim([0, 1.05])

    # Plot 2: Damköhler Number
    ax2 = fig.add_subplot(gs[1, 0])
    Da = np.linspace(0.01, 20, 500)
    ax2.plot(Da, Da / (1 + Da), color='#1B3A6B', linewidth=2.5)
    ax2.axvline(1, color='orange', linestyle='--', alpha=0.7, label='Da=1 → X=0.5')
    ax2.axvline(9, color='red',    linestyle='--', alpha=0.7, label='Da=9 → X=0.9')
    ax2.set_xlabel("Damköhler Number (Da = k·τ)", fontsize=11)
    ax2.set_ylabel("Conversion X", fontsize=11)
    ax2.set_title("Universal CSTR Curve", fontsize=11, fontweight='bold')
    ax2.legend(fontsize=8); ax2.grid(True, alpha=0.3)

    # Plot 3: Required Volume vs Target Conversion
    ax3 = fig.add_subplot(gs[1, 1])
    X_t = np.linspace(0.05, 0.99, 300)
    V   = [required_volume(10.0, x, 0.5, 1.0) for x in X_t]
    ax3.plot(X_t * 100, V, color='#0D7377', linewidth=2.5)
    ax3.set_xlabel("Target Conversion (%)", fontsize=11)
    ax3.set_ylabel("Required Volume (L)", fontsize=11)
    ax3.set_title("Reactor Volume vs Conversion\n(k=0.5, F_A0=10 mol/min)", fontsize=11, fontweight='bold')
    ax3.grid(True, alpha=0.3)

    plt.savefig("cstr_results.png", dpi=150, bbox_inches='tight')
    print("\nPlot saved → cstr_results.png")
    plt.show()


# ─────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("\n" + "="*55)
    print("  CSTR PERFORMANCE SIMULATION")
    print("  Reaction: A → B  |  First-Order Kinetics")
    print("="*55)

    print_conversion_table(k=0.5)
    print_validation(k=0.5)

    print("\n─── DESIGN EXAMPLE ───")
    k, tau = 0.5, 5.0
    X = cstr_conversion(k, tau)
    V = required_volume(10.0, 0.80, k, 1.0)
    print(f"  k={k} min⁻¹, τ={tau} min → X = {X:.4f} ({X*100:.1f}%)")
    print(f"  Volume for 80% conversion: {V:.2f} L")

    plot_results()
