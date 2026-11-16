# -*- coding: utf-8 -*-
"""
PHASE VIII -- QAOA EMULATION ON THE FOUR-OPTION PLANT PROBLEM WITH TWO LIMITS (CAPITAL AND SCRAP)

How to run in Google Colab
  1. Run phase2_route_estimates_colab.py, phase3_route_model_colab.py, phase6_mass_energy_balance_colab.py
     and phase7_two_limit_selection_colab.py first (same runtime), or upload the folder outputs_phase3
     (needs plant_inputs.csv, scenario_table.csv, mass_energy_balance.csv) and, if available, outputs_phase4
     (twolimit_plant_options.csv and twolimit_decisions.csv from phase 7).
  2. Run:  !python -X utf8 phase8_twolimit_qaoa_colab.py
  Needs only numpy and scipy (no Qiskit). Run time: see the last printed line.
  If outputs_phase4/twolimit_decisions.csv (from the exact classical solver, phase 7) is present, the best
  plans found here by brute force must equal that file in EVERY case (assert); if it is absent, a note is
  printed. If outputs_phase4/twolimit_plant_options.csv is present, the scrap mass per tonne used here must
  equal the values in that file (assert).
  Scrap mass per tonne of steel: read from the iron balance of phase 6 (stream "scrap (iron-balance scaling
  M)" of outputs_phase3/mass_energy_balance.csv) for s = 0, 0.50, 0.95, with the same reading code as phase 7.

Problem (one plan = one option per plant)
  Four plants (Jamshedpur, Vijayanagar, Dolvi, Bhilai). Each plant takes ONE option:
  N = do not convert, S0 / S50 / S95 = convert to hydrogen direct reduction + electric arc furnace with a scrap
  share of 0, 0.50 or 0.95. Maximise the sum of abatement, subject to a capital limit B and a scrap limit S.
  Encoding: 2 bits per plant (00 = N, 01 = S0, 10 = S50, 11 = S95), plant J first -> 8 bits, 256 plans.
  Every 8-bit string is a valid plan, so no "one option per plant" penalty is needed.

Methods (all are exact state-vector emulations with numpy on a classical computer)
  G       Grover-mixer QAOA. Start = equal superposition |F> of the plans that satisfy BOTH limits.
          Phase separator exp(-i gamma H_n), H_n = minus abatement rescaled to 0..1 over the feasible plans.
          Mixer 1 - (1 - exp(-i beta)) |F><F|. The state never leaves F (asserted, 1e-9).
  Three PENALTY methods: standard QAOA (X mixer on all 8 qubits, start = equal superposition of all 256
  plans) on an energy H = minus abatement + penalty, rescaled to 0..1. Capital and scrap are divided by their
  own limit, so the violations have no units and the two limits are comparable; one penalty weight P for both.
    quadratic form:  H = minus abatement + P x [max(0, capital/B - 1)^2 + max(0, scrap/S - 1)^2]
    hinge form:      H = minus abatement + P x [max(0, capital/B - 1)   + max(0, scrap/S - 1)]   (linear)
    [ASSUMPTION: two common forms; results for penalty methods depend on this choice]
  O_fixed  quadratic form with a FIXED PENALTY RULE, P = 1 + sum over plants of the largest abatement of that
           plant. This rule DOES NOT ALWAYS SEPARATE: with it, an infeasible plan can have a lower energy than
           the best feasible plan. The script CHECKS this in every case and writes the answer in column
           penalty_separates (it is NOT assumed to hold).
  O_min    quadratic form, P = 1.05 x the smallest P for which every infeasible plan has a higher energy than
           the best feasible plan (computed exactly from the 256 plans).
  O_hinge  hinge form, P = 1.05 x the smallest P for which every infeasible plan has a higher energy than the
           best feasible plan (computed exactly from the 256 plans).
           O_min and O_hinge are marked "not needed" (and not run) in a case where no infeasible plan has more
           abatement than the best feasible plan, because then the smallest such P is zero.
  In the summary, the means of each penalty method are given twice: over the cases where its penalty separates
  (with the number of such cases), and over all cases in which the method was run.
  Reference lines: equal chance over the feasible plans (1/|F|) and over all 256 plans (1/256).
  A slack-variable QUBO is NOT built here: two inequality limits with non-integer coefficients would need
  many extra slack qubits on top of the 8 decision qubits. This script makes no statement about that encoding.

Reported for each method and depth
  P_best        probability of the best plan (one plan; ties in abatement are broken by lower capital, then
                lower scrap; the number of feasible plans tied on abatement is written to n_tied_best)
  TWO COUNTS are given in the summary, and both must be quoted together:
                n_most_likely_is_best            the single best plan (tie rule) is the most likely output
                n_most_likely_has_best_abatement a plan with the best abatement is the most likely output
                The second is the headline count, because a plan tied on abatement is as good as the best plan.
  P_feasible    probability of a plan that satisfies both limits
  most likely plan, and whether it is the best plan
  P_best_abatement  probability of ANY feasible plan whose abatement equals the best abatement. With renewable
                electricity the route intensity does not depend on the scrap share, so S0, S50 and S95 of a plant
                give the same abatement and many plans tie. The Grover mixer gives tied plans the same
                probability, so there "most likely plan = best plan" depends on the order of the plans; use
                most_likely_has_best_abatement (1 if no plan is more likely than the tied best plans) instead.
  exp_ratio     expected abatement divided by the best abatement. RULE: an infeasible output counts as ZERO
                abatement (it cannot be built). A feasible output with negative abatement counts as negative.

What this script does NOT show (read before citing any number)
  * The feasible plans are listed on a classical computer BEFORE the emulation (all 256 plans are evaluated).
    That listing already solves a problem of this size. The best plan is known before any QAOA angle is set.
    The penalty weights of O_min and O_hinge are also computed from that full listing and from the known best plan.
  * No circuit is built for preparing |F>, for the Grover mixer, or for the penalty arithmetic
    max(0, .) and max(0, .)^2. Their gate cost is not estimated.
  * QAOA angles are optimised on the exact state vector, with the exact energy, which a real device cannot do.
  * Nothing here shows a speed-up or any advantage over the classical brute force used in this same script.
  * Depths: p = 1..6 on the main case, p = 1..3 on the sweep.

Tags: [VERIFIED] read in the source document (PDF page and quote);
      [ASSUMPTION] study choice.

Outputs (outputs_phase4/)
  twolimit_qaoa_main.csv      main case (India grid, f_B = 0.50, f_S = 0.25), method x depth 1..6
  twolimit_qaoa_sweep.csv     both electricity cases x 10 capital limits x 5 scrap limits, depth 1..3;
                              cases without a real choice are listed with method = not_run
  twolimit_qaoa_summary.csv   per method and depth: counts, mean and minimum P_best, and so on; for the penalty
                              methods once over the cases where the penalty separates and once over all cases run
  software_versions_phase8_twolimit_qaoa_colab.csv   python, numpy and scipy versions of the run
Tie rule and rounding: the same as phase 7 (abatement, capital and scrap rounded to 9 decimals; highest
abatement, then lower capital, then lower scrap).
"""
import os, csv, math, time, platform
import numpy as np
import scipy
from scipy.optimize import minimize

T0 = time.time()
OUT = "outputs_phase4"
os.makedirs(OUT, exist_ok=True)
SEED = 12345
print("python", platform.python_version(), "| numpy", np.__version__, "| scipy", scipy.__version__)

def z0(v, nd=6):   # format without printing "-0.000000"
    v = round(float(v), nd)
    return f"{(0.0 if v == 0 else v):.{nd}f}"

# =============================================================================
# PART 1 -- CONSTANTS
# =============================================================================
# Vogl, Ahman & Nilsson (2018), J. Clean. Prod. 203:736-745; file 1_s2_0_S0959652618326301_main_pdf.txt
SEC_H2DRI_0SCRAP = 3.48    # MWh/t liquid steel. PDF PAGE 4: "A specific energy consumption (SEC) of 3.48 MWh/t of
                           # liquid"                                                              [VERIFIED]
SEC_EAF_100SCRAP = 0.667   # MWh/t. PDF PAGE 5: "less energy (0.667 MWh/tLS) than from pure DRI"  [VERIFIED]
DIRECT_CO2       = 0.053   # t CO2/t. PDF PAGE 6: "graphite electrodes would result in emissions of 53 kg CO2
                           # per"                                                                 [VERIFIED]
REST_FUEL_GJ     = 2.2     # GJ/t. PDF PAGE 4: "For pelletizing and rolling an additional 2.2 GJ/t of fuel as"
REST_ELEC_MWH    = 0.140   # MWh/t. PDF PAGE 5: "well as 140 kWh/t electricity are required"      [VERIFIED both]
CAPEX_ELECTROLYSER = 160.0 # EUR per t of yearly capacity. PDF PAGE 6: "It includes 160 EUR/t capacity for the
CAPEX_SHAFT        = 230.0 # electrolyser, 230 EUR/t for the shaft and 184 EUR/t for the EAF"     [VERIFIED]
CAPEX_EAF          = 184.0 #                                                                      [VERIFIED]
CAPEX_TOTAL        = 574.0 # PDF PAGE 6: "The CAPEX for a H-DR plant was calculated as 574 EUR per tonne
                           # capacity"                                                            [VERIFIED]
assert CAPEX_ELECTROLYSER + CAPEX_SHAFT + CAPEX_EAF == CAPEX_TOTAL
SCRAP_KG_PER_T_AT_50 = 536.0 # kg scrap per t steel at 50% scrap. PDF PAGE 4: "equal shares of HBI and scrap,
                             # 738 kg of pellets, 536 kg of scrap and 25 kg of hydrogen are needed per tonne
                             # of steel"                                                          [VERIFIED]
# The 536 kg is used ONLY as a check. The scrap mass per tonne of steel is READ from the iron balance of
# phase6_mass_energy_balance_colab.py (below). That balance uses the compositions printed by Vogl et al.
# (PDF PAGE 2: "Scrap charged to the EAF contains 95% iron and 5% inert substances."; PDF PAGE 3: "The liquid
# steel product contains only iron") and closes the iron balance at 1000 kg iron per tonne of steel at every
# scrap share (by construction of that scaling). The second scaling of phase 6 ("proportional scaling P", 1072 kg x s)
# does not close the iron balance away from s = 0.50 and is NOT used.
# In phase 6 the scrap share s is read as scrap mass / (scrap mass + HBI mass) of the furnace charge. [ASSUMPTION]
SCRAP_STREAM = "scrap (iron-balance scaling M)"          # stream name as written by phase 6
SCRAP_IRON_STREAM = "iron in scrap (M)"
X_FE_SCRAP = 0.95            # iron share of scrap, Vogl PDF PAGE 2 (quote above)                 [VERIFIED]
BAL_CSV = "outputs_phase3/mass_energy_balance.csv"
# Electrolyser and shaft capacity scale with the ore share (1 - s); plant output = capacity.     [ASSUMPTION]
# SEC is linear in s between Vogl's two end points (as in phase3_route_model_colab.py).          [ASSUMPTION]

# IPCC 2006 Guidelines, Vol. 2, Ch. 1, Table 1.4; PDF PAGE 24:
# "Natural Gas 15.3 1 56 100 54 300 58 300" -> 56,100 kg CO2/TJ                                   [VERIFIED]
EF_NATGAS = 0.0561         # t CO2/GJ
                           # the 2.2 GJ fuel is taken as natural gas; Vogl does not name the fuel [ASSUMPTION]

# CEA, CO2 Baseline Database for the Indian Power Sector, User Guide Version 20.0;
# PDF PAGE 24: "2023-24* 1204.51 1408.45 225.83 1655.70 0.727"
# and PDF PAGE 16: "Including cross-border power transfers 0.727 0.962 0.552 0.757"               [VERIFIED]
GRID_INDIA     = 0.727     # t CO2/MWh
GRID_RENEWABLE = 0.0       # t CO2/MWh, operation only                                            [ASSUMPTION]
ELECTRICITY = [("India grid", GRID_INDIA), ("Renewable", GRID_RENEWABLE)]

# EY-Parthenon, WWF-India & CII-GBC (February 2026); file ey_closing_the_loop_scrap_markets_to_power_indias_
# green_stee.txt, PDF PAGE 18: "In 2024, India's total ferrous scrap consumption stood at 34.2 million tons."
S_NATIONAL = 34.2          # Mt scrap per year (2024)                                             [VERIFIED]

OPTIONS      = ["N", "S0", "S50", "S95"]            # option code 0, 1, 2, 3 = bits 00, 01, 10, 11
OPTION_CHAR  = ["N", "0", "5", "9"]
OPTION_SCRAP = [None, 0.0, 0.50, 0.95]              # scrap shares of the study                   [ASSUMPTION]
F_B_GRID = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]                         # [ASSUMPTION]
F_S_GRID = [0.10, 0.25, 0.50, 0.75, 1.00]                                                       # [ASSUMPTION]
MAIN_CASE = ("India grid", 0.50, 0.25)                                                          # [ASSUMPTION]
P_MAX_MAIN, P_MAX_SWEEP = 6, 3
N_RANDOM = 4               # random starts per depth                                              [ASSUMPTION]

def sec(s):  return SEC_H2DRI_0SCRAP - s * (SEC_H2DRI_0SCRAP - SEC_EAF_100SCRAP)
def I_new(s, g): return DIRECT_CO2 + sec(s) * g + REST_FUEL_GJ * EF_NATGAS + REST_ELEC_MWH * g
def capex(s): return CAPEX_EAF + (1.0 - s) * (CAPEX_ELECTROLYSER + CAPEX_SHAFT)


def read_scrap_t_per_t(path):
    """Scrap mass per tonne of steel (t/t) for s = 0, 0.50, 0.95, read from the iron balance of phase 6."""
    with open(path, encoding="utf-8") as f:
        bal = list(csv.DictReader(f))
    scr = {round(float(r["scrap_share"]), 2): float(r["value"]) for r in bal if r["stream"] == SCRAP_STREAM}
    fe = {round(float(r["scrap_share"]), 2): float(r["value"]) for r in bal if r["stream"] == SCRAP_IRON_STREAM}
    for s in (0.0, 0.50, 0.95):
        if s not in scr or s not in fe:
            raise KeyError(f"{path} has no row '{SCRAP_STREAM}' / '{SCRAP_IRON_STREAM}' for scrap share {s}")
        assert scr[s] >= 0 and abs(fe[s] - X_FE_SCRAP * scr[s]) < 1e-3, "scrap row is not 95 % iron"
        assert fe[s] <= 1000.0 + 1e-6, f"iron in scrap exceeds the iron in one tonne of steel at s = {s}"
    assert scr[0.0] == 0.0 and scr[0.0] < scr[0.50] < scr[0.95]
    assert abs(scr[0.50] / SCRAP_KG_PER_T_AT_50 - 1) <= 0.01, \
        "scrap at 50 % scrap from the iron balance differs from Vogl's printed 536 kg by more than 1 %"
    return {s: scr[s] / 1000.0 for s in (0.0, 0.50, 0.95)}


if not os.path.exists(BAL_CSV):
    raise FileNotFoundError("Run phase6_mass_energy_balance_colab.py first, or upload the file: " + BAL_CSV)
SCRAP_T_PER_T = read_scrap_t_per_t(BAL_CSV)
def m_scrap(s): return SCRAP_T_PER_T[round(s, 2)]
assert abs(capex(0.0) - CAPEX_TOTAL) < 1e-12
print(f"Scrap t per t steel (iron balance, phase 6): s = 0 -> {m_scrap(0.0):.6f}, 0.50 -> {m_scrap(0.5):.6f}, "
      f"0.95 -> {m_scrap(0.95):.6f}")

# the scrap values must be the ones used by the exact classical solver (phase 7), when its table is present
OPT_CSV = f"{OUT}/twolimit_plant_options.csv"
if os.path.exists(OPT_CSV):
    with open(OPT_CSV, encoding="utf-8") as f:
        _n = 0
        for r in csv.DictReader(f):
            if r["option"] in ("S0", "S50", "S95"):
                assert abs(float(r["scrap_t_per_t"]) - m_scrap(float(r["scrap_share"]))) < 5e-7, \
                    f"scrap mass per tonne differs from {OPT_CSV} ({r['plant']}, {r['option']})"
                _n += 1
    assert _n == 24, f"{OPT_CSV}: expected 24 converted-option rows, found {_n}"
    print(f"  equal to the scrap values in {OPT_CSV} ({_n} rows checked)")
else:
    print(f"NOTE: {OPT_CSV} not found -> scrap values were NOT compared with the phase 7 table.")

# cross-check of the route model against the verified Phase III output (not a self-check: read from the CSV)
with open("outputs_phase3/scenario_table.csv", encoding="utf-8") as f:
    P3 = {(float(r["scrap_share"]), r["electricity"]): float(r["I_new_tCO2_per_t"]) for r in csv.DictReader(f)}
for s_chk in (0.95, 0.50, 0.0):
    for e_chk, g_chk in ELECTRICITY:
        # the CSV holds 6 decimals; I(0.50, 0.727) = 1.7856345 lies on a rounding edge, hence 5e-7 + 1e-9
        assert abs(I_new(s_chk, g_chk) - P3[(s_chk, e_chk)]) <= 5e-7 + 1e-9, \
            f"I({s_chk}, {e_chk}) differs from outputs_phase3/scenario_table.csv"
print(f"Route model check vs Phase III: I(0.95, 0.727) = {I_new(0.95, GRID_INDIA):.6f}, "
      f"I(0, 0.727) = {I_new(0.0, GRID_INDIA):.6f}  (match)")

# plants: output and GEI of the 2026 revised draft, read from the plant table written by phase 3
with open("outputs_phase3/plant_inputs.csv", encoding="utf-8") as f:
    PLANTS = list(csv.DictReader(f))
SHORT = ["J", "V", "D", "B"]
EXPECT = ["Jamshedpur", "Vijayanagar", "Dolvi", "Bhilai"]
assert len(PLANTS) == 4 and all(e in r["plant"] for e, r in zip(EXPECT, PLANTS)), "plant order must be J, V, D, B"
Q   = np.array([float(r["output_2023_24_t_rev2026"]) for r in PLANTS])
GEI = np.array([float(r["GEI_2023_24_rev2026"]) for r in PLANTS])
NP_, NQ, NPLAN = 4, 8, 256
C_ALL0 = float(Q.sum() * CAPEX_TOTAL / 1e6)          # million EUR: all four plants converted at s = 0
print(f"C_all0 = {C_ALL0:.3f} million EUR; S_nat = {S_NATIONAL} Mt/yr")

# =============================================================================
# PART 2 -- TABLES OVER ALL 256 PLANS AND BRUTE-FORCE OPTIMUM
# =============================================================================
# plan index k (0..255): the 8-bit binary string of k, read from the left, is plant J first.
OPT_OF = np.array([[(k >> (2 * (NP_ - 1 - i))) & 3 for i in range(NP_)] for k in range(NPLAN)])   # 256 x 4
def plan_label(k): return "".join(OPTION_CHAR[o] for o in OPT_OF[k])
assert plan_label(0) == "NNNN" and plan_label(0b11000010) == "9NN5" and len({plan_label(k) for k in range(NPLAN)}) == 256

def option_tables(g):
    """per plant and option: abatement (Mt CO2e/yr), capital (million EUR), scrap (Mt/yr)"""
    A = np.zeros((NP_, 4)); C = np.zeros((NP_, 4)); U = np.zeros((NP_, 4))
    for o in (1, 2, 3):
        s = OPTION_SCRAP[o]
        A[:, o] = Q * (GEI - I_new(s, g)) / 1e6
        C[:, o] = Q * capex(s) / 1e6
        U[:, o] = Q * m_scrap(s) / 1e6
    return A, C, U

def plan_tables(g):
    A, C, U = option_tables(g); idx = np.arange(NP_)[None, :]
    return A[idx, OPT_OF].sum(1), C[idx, OPT_OF].sum(1), U[idx, OPT_OF].sum(1), A

TOL = 1e-9        # same tolerance as phase 7 (limits and ties)

def best_plan(abv, capv, scrv, feas):
    """Tie rule and rounding of phase 7: highest abatement, then lower capital, then lower scrap, each rounded to
    9 decimals. Returns the plan index and the number of feasible plans tied with it on abatement (within TOL)."""
    f = np.where(feas)[0]
    order = np.lexsort((np.round(scrv[f], 9), np.round(capv[f], 9), -np.round(abv[f], 9)))
    k = int(f[order[0]])
    return k, int(np.sum(np.abs(abv[f] - abv[k]) <= TOL))

# =============================================================================
# PART 3 -- STATE-VECTOR QAOA
# =============================================================================
def fwht(a):
    h, n = 1, a.size
    while h < n:
        a = a.reshape(-1, 2, h)
        a = np.stack([a[:, 0] + a[:, 1], a[:, 0] - a[:, 1]], axis=1).reshape(-1)
        h *= 2
    return a / np.sqrt(n)

def zsum(m):
    return np.array([sum(1 - 2 * ((k >> q) & 1) for q in range(m)) for k in range(2 ** m)], dtype=float)

def rescale(d):
    lo, hi = d.min(), d.max()
    return (d - lo) / (hi - lo) if hi > lo else np.zeros_like(d)

_t = np.random.default_rng(1).normal(size=NPLAN)
assert np.allclose(fwht(fwht(_t)), _t, atol=1e-12), "the Walsh-Hadamard transform must be its own inverse"
Z8 = zsum(NQ)

def qaoa_state(params, p, Hn, Zs):
    g, b = params[:p], params[p:]
    psi = np.full(Hn.size, 1 / np.sqrt(Hn.size), dtype=complex)
    for l in range(p):
        psi = psi * np.exp(-1j * g[l] * Hn)
        psi = fwht(np.exp(-1j * b[l] * Zs) * fwht(psi))          # exp(-i*beta*sum X) in the Hadamard basis
    return psi

def qaoa_energy(params, p, Hn, Zs):
    psi = qaoa_state(params, p, Hn, Zs)
    return float(np.real(np.vdot(psi, Hn * psi)))

def _starts(p, prev, rng, hi):
    s = [np.concatenate([np.linspace(0.1, 1.0, p) * 3.0, np.linspace(1.0, 0.1, p) * 0.75])]   # linear ramp
    if prev is not None and len(prev) == 2 * (p - 1):                                         # previous depth
        g, b = prev[:p - 1], prev[p - 1:]
        s.append(np.concatenate([g, [g[-1]], b, [b[-1]]]))
    return s + [rng.uniform(0, hi, 2 * p) for _ in range(N_RANDOM)]                           # random starts

def qaoa_x_mixer(H, pmax, seed):
    """standard QAOA (X mixer) on a diagonal energy over the 256 plans; probability vector for p = 1..pmax"""
    Hn = rescale(H); rng_l = np.random.default_rng(seed)
    out, prev = [], None
    for p in range(1, pmax + 1):
        best = min((minimize(qaoa_energy, s0, args=(p, Hn, Z8), method="L-BFGS-B")
                    for s0 in _starts(p, prev, rng_l, np.pi)), key=lambda r: r.fun)
        prev = best.x; pr = np.abs(qaoa_state(best.x, p, Hn, Z8)) ** 2
        assert abs(pr.sum() - 1) < 1e-9
        out.append(pr)
    return out

def grover_state(prm, p, fn):
    m = fn.size; s = np.full(m, 1 / np.sqrt(m), dtype=complex); psi = s.copy()
    for l in range(p):
        psi = psi * np.exp(-1j * prm[l] * fn)
        psi = psi + (np.exp(-1j * prm[p + l]) - 1) * s * np.vdot(s, psi)     # 1 - (1 - exp(-i*beta)) |F><F|
    return psi

def grover_energy(prm, p, fn):
    psi = grover_state(prm, p, fn)
    return float(np.real(np.vdot(psi, fn * psi)))

def qaoa_grover_mixer(abv, feas, pmax, seed):
    """Grover-mixer QAOA on the feasible plans; probability over ALL 256 plans (0 outside F)"""
    f = np.where(feas)[0]; fn = rescale(-abv[f]); rng_l = np.random.default_rng(seed)
    out, prev = [], None
    for p in range(1, pmax + 1):
        best = min((minimize(grover_energy, s0, args=(p, fn), method="L-BFGS-B")
                    for s0 in _starts(p, prev, rng_l, 2 * np.pi)), key=lambda r: r.fun)
        prev = best.x; pr = np.zeros(abv.size); pr[f] = np.abs(grover_state(best.x, p, fn)) ** 2
        assert abs(pr[f].sum() - 1.0) < 1e-9 and pr[~feas].sum() == 0.0, "Grover-mixer state left the feasible set"
        out.append(pr)
    return out

def outcome(prob, abv, feas, k_best):
    """P(best plan), P(feasible), most likely plan, expected abatement / best (infeasible outputs count as zero)"""
    tie = feas & (np.abs(abv - abv[k_best]) <= TOL)          # feasible plans with the same abatement as the best plan
    k_ml = int(np.argmax(prob))
    return (float(prob[k_best]), float(prob[feas].sum()), k_ml, float(prob[feas] @ abv[feas] / abv[k_best]),
            float(prob[tie].sum()), int(prob[k_ml] <= prob[tie].max() + 1e-12))

# =============================================================================
# PART 4 -- ONE CASE
# =============================================================================
def setup_case(elec, g, fB, fS):
    abv, capv, scrv, A_opt = plan_tables(g)
    B, S = fB * C_ALL0, fS * S_NATIONAL
    feas = (capv <= B + TOL) & (scrv <= S + TOL)
    assert feas[0], "the plan NNNN must always be feasible"
    k_best, n_tied = best_plan(abv, capv, scrv, feas)
    # Penalty forms            [ASSUMPTION: two common forms; results for penalty methods depend on this choice]
    over_c, over_s = np.maximum(0.0, capv / B - 1.0), np.maximum(0.0, scrv / S - 1.0)   # normalised violations
    viol = over_c ** 2 + over_s ** 2          # quadratic form (O_fixed, O_min)
    viol_h = over_c + over_s                  # hinge form, linear in the normalised violation (O_hinge)
    for v_ in (viol, viol_h):
        if not (np.all(v_[feas] == 0) and np.all(v_[~feas] > 0)):
            raise ValueError("a penalty term is not zero on a feasible plan or not positive on an infeasible plan")
    P_fixed = 1.0 + float(A_opt.max(axis=1).sum())          # fixed penalty rule; does not always separate
    inf = ~feas
    if inf.any():
        # smallest P for which no infeasible plan has a LOWER energy than the best feasible plan (equal at P itself)
        P_need = float(np.max((abv[inf] - abv[k_best]) / viol[inf]))
        P_need_h = float(np.max((abv[inf] - abv[k_best]) / viol_h[inf]))
        fixed_sep = bool(np.all(-abv[inf] + P_fixed * viol[inf] > -abv[k_best]))
        min_viol = float(np.sqrt(viol[inf].min()))
    else:
        P_need, P_need_h, fixed_sep, min_viol = 0.0, 0.0, True, float("nan")
    # O_min / O_hinge are "not needed" when no infeasible plan has more abatement than the best feasible plan
    P_min = 1.05 * P_need if P_need > 0 else None
    P_hinge = 1.05 * P_need_h if P_need_h > 0 else None
    if P_min is not None:
        assert np.all(-abv[inf] + P_min * viol[inf] > -abv[k_best]), "O_min penalty does not separate"
        assert np.all(-abv[inf] + P_hinge * viol_h[inf] > -abv[k_best]), "O_hinge penalty does not separate"
    return dict(elec=elec, g=g, fB=fB, fS=fS, B=B, S=S, abv=abv, capv=capv, scrv=scrv, feas=feas, viol=viol,
                viol_h=viol_h, k_best=k_best, n_tied=n_tied, nF=int(feas.sum()), P_fixed=P_fixed, P_need=P_need,
                P_need_h=P_need_h, P_min=P_min, P_hinge=P_hinge, fixed_sep=fixed_sep, min_viol=min_viol,
                real_choice=bool(k_best != 0 and feas.sum() >= 2))

def run_case(c, pmax):
    """returns rows: one per method and depth"""
    rows = []
    def add(method, P, sep, probs):
        print(f"    done {c['elec']} f_B={c['fB']:.2f} f_S={c['fS']:.2f} {method}  t = {time.time() - T0:.0f} s", flush=True)
        for p, pr in enumerate(probs, 1):
            pb, pf, k_ml, ratio, pt, ml_tie = outcome(pr, c["abv"], c["feas"], c["k_best"])
            rows.append(dict(method=method, penalty=P, penalty_separates=sep, p=p, P_best=pb, P_feasible=pf,
                             most_likely_plan=plan_label(k_ml), most_likely_is_best=int(k_ml == c["k_best"]),
                             exp_ratio=ratio, P_tied=pt, ml_tied=ml_tie))
    add("G", None, "", qaoa_grover_mixer(c["abv"], c["feas"], pmax, SEED))
    add("O_fixed", c["P_fixed"], int(c["fixed_sep"]), qaoa_x_mixer(-c["abv"] + c["P_fixed"] * c["viol"], pmax, SEED))
    if c["P_min"] is not None:
        add("O_min", c["P_min"], 1, qaoa_x_mixer(-c["abv"] + c["P_min"] * c["viol"], pmax, SEED))
        add("O_hinge", c["P_hinge"], 1, qaoa_x_mixer(-c["abv"] + c["P_hinge"] * c["viol_h"], pmax, SEED))
    return rows

# =============================================================================
# PART 5 -- MAIN CASE (depth 1..6)
# =============================================================================
g_main = dict(ELECTRICITY)[MAIN_CASE[0]]
M = setup_case(MAIN_CASE[0], g_main, MAIN_CASE[1], MAIN_CASE[2])
kb = M["k_best"]
print(f"\nMAIN CASE: {M['elec']}, f_B = {M['fB']}, f_S = {M['fS']}: B = {M['B']:.3f} million EUR, S = {M['S']:.3f} Mt/yr")
print(f"  feasible plans |F| = {M['nF']} of 256; best plan {plan_label(kb)}: abatement {M['abv'][kb]:.6f} Mt CO2e/yr, "
      f"capital {M['capv'][kb]:.3f}, scrap {M['scrv'][kb]:.6f}; plans tied on abatement: {M['n_tied']}")
print(f"  fixed penalty rule (O_fixed, quadratic form) P = {M['P_fixed']:.6f}; separates every infeasible plan from the "
      f"best plan: {M['fixed_sep']}")
print(f"  smallest violation that occurs (root of the squared sum, as a share of the limit) = {M['min_viol']:.6f}")
print(f"  smallest separating penalty, quadratic form = {M['P_need']:.6f}; O_min uses 1.05 x = "
      f"{'not needed' if M['P_min'] is None else format(M['P_min'], '.6f')}")
print(f"  smallest separating penalty, hinge form = {M['P_need_h']:.6f}; O_hinge uses 1.05 x = "
      f"{'not needed' if M['P_hinge'] is None else format(M['P_hinge'], '.6f')}")
assert M["real_choice"], "main case has no real choice"
t1 = time.time()
MAIN_ROWS = run_case(M, P_MAX_MAIN)
print(f"  (main case run time {time.time() - t1:.1f} s)")
EQ_F, EQ_ALL = 1.0 / M["nF"], 1.0 / NPLAN
EQ_RATIO_F = float(M["abv"][M["feas"]].mean() / M["abv"][kb])            # expected ratio of an equal-chance draw from F
EQ_RATIO_ALL = float(M["abv"][M["feas"]].sum() / NPLAN / M["abv"][kb])   # ... from all 256 (infeasible = zero)
HEAD_M = ["electricity", "f_B", "f_S", "B_million_EUR", "S_Mt_scrap", "n_feasible", "best_plan", "best_abatement_Mt",
          "method", "penalty_P1_P2", "penalty_separates", "p", "P_best", "P_feasible", "most_likely_plan",
          "most_likely_is_best", "exp_abatement_ratio", "equal_chance_feasible", "equal_chance_all_256",
          "n_tied_best", "P_best_abatement", "most_likely_has_best_abatement"]
def fmt_rows(c, rows):
    for r in rows:
        yield [c["elec"], f"{c['fB']:.2f}", f"{c['fS']:.2f}", z0(c["B"]), z0(c["S"]), c["nF"], plan_label(c["k_best"]),
               z0(c["abv"][c["k_best"]]), r["method"], "" if r["penalty"] is None else z0(r["penalty"]),
               r["penalty_separates"], r["p"], z0(r["P_best"]), z0(r["P_feasible"]), r["most_likely_plan"],
               r["most_likely_is_best"], z0(r["exp_ratio"]), z0(1.0 / c["nF"]), z0(1.0 / NPLAN),
               c["n_tied"], z0(r["P_tied"]), r["ml_tied"]]
with open(f"{OUT}/twolimit_qaoa_main.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f); wr.writerow(HEAD_M)
    for row in fmt_rows(M, MAIN_ROWS):
        wr.writerow(row)
    for name, pb, pf, ratio in [("equal_chance_feasible", EQ_F, 1.0, EQ_RATIO_F),
                                ("equal_chance_all_256", EQ_ALL, M["nF"] / NPLAN, EQ_RATIO_ALL)]:
        wr.writerow([M["elec"], f"{M['fB']:.2f}", f"{M['fS']:.2f}", z0(M["B"]), z0(M["S"]), M["nF"], plan_label(kb),
                     z0(M["abv"][kb]), name, "", "", "", z0(pb), z0(pf), "", "", z0(ratio), z0(EQ_F), z0(EQ_ALL),
                     M["n_tied"], z0(pb * M["n_tied"]), ""])
print(f"  reference: equal chance over F = {EQ_F:.6f} (expected ratio {EQ_RATIO_F:.6f}); over all 256 = {EQ_ALL:.6f} "
      f"(P feasible {M['nF'] / NPLAN:.6f}, expected ratio {EQ_RATIO_ALL:.6f})")
print("  method  p   P_best  P_feasible  most likely  is best  exp. ratio")
for r in MAIN_ROWS:
    print(f"  {r['method']:7s} {r['p']}  {r['P_best']:.6f}  {r['P_feasible']:.6f}   {r['most_likely_plan']}       "
          f"{r['most_likely_is_best']}     {r['exp_ratio']:.6f}")

# =============================================================================
# PART 6 -- SWEEP (depth 1..3): both electricity cases x all f_B x all f_S
# =============================================================================
print(f"\nSWEEP: {len(ELECTRICITY)} x {len(F_B_GRID)} x {len(F_S_GRID)} = "
      f"{len(ELECTRICITY) * len(F_B_GRID) * len(F_S_GRID)} cases, depth 1..{P_MAX_SWEEP}")
CASES, t1 = [], time.time()
seen_literal, seen_F = {}, {}
for elec, g in ELECTRICITY:
    for fB in F_B_GRID:
        for fS in F_S_GRID:
            c = setup_case(elec, g, fB, fS)
            c["id"] = f"{elec}|{fB:.2f}|{fS:.2f}"
            # (a) duplicate as defined in the task: same abatement, capital and scrap tables AND same limits
            key_a = (c["abv"].round(9).tobytes(), c["capv"].round(9).tobytes(), c["scrv"].round(9).tobytes(),
                     round(c["B"], 9), round(c["S"], 9))
            c["dup_literal_of"] = seen_literal.setdefault(key_a, c["id"])
            # (b) same abatement table and same feasible set: same best plan and the SAME Grover-mixer problem;
            #     the penalty energies still differ with the limits
            key_b = (c["abv"].round(9).tobytes(), c["feas"].tobytes())
            c["same_F_as"] = seen_F.setdefault(key_b, c["id"])
            c["rows"] = run_case(c, P_MAX_SWEEP) if c["real_choice"] else []
            CASES.append(c)
print(f"  (sweep run time {time.time() - t1:.1f} s)")
n_real = sum(c["real_choice"] for c in CASES)
n_lit = sum(c["real_choice"] and c["dup_literal_of"] == c["id"] for c in CASES)
n_dF = sum(c["real_choice"] and c["same_F_as"] == c["id"] for c in CASES)
print(f"  cases with a real choice: {n_real} of {len(CASES)}; distinct by tables + limits: {n_lit}; "
      f"distinct by abatement table + feasible set: {n_dF}")
print(f"  fixed penalty rule (O_fixed) separates in {sum(c['fixed_sep'] for c in CASES if c['real_choice'])} of {n_real} "
      f"real-choice cases; O_min and O_hinge run in {sum(c['P_min'] is not None for c in CASES if c['real_choice'])} "
      f"(not needed in the others)")

# the main case inside the sweep must reproduce depths 1..3 of PART 5 (same seed, same sequence of starts)
cm = next(c for c in CASES if (c["elec"], c["fB"], c["fS"]) == MAIN_CASE)
for r in cm["rows"]:
    m = next(x for x in MAIN_ROWS if x["method"] == r["method"] and x["p"] == r["p"])
    assert abs(m["P_best"] - r["P_best"]) < 1e-12, "sweep and main case disagree"

HEAD_S = ["electricity", "f_B", "f_S", "B_million_EUR", "S_Mt_scrap", "n_feasible", "best_plan", "best_abatement_Mt",
          "n_tied_best", "real_choice", "distinct_problem", "same_feasible_set_as", "distinct_feasible_set",
          "fixed_penalty", "fixed_penalty_separates", "smallest_separating_penalty_quadratic",
          "smallest_separating_penalty_hinge", "O_min_and_O_hinge_status", "smallest_violation_share",
          "method", "penalty_P1_P2", "p", "P_best", "P_feasible", "most_likely_plan", "most_likely_is_best",
          "exp_abatement_ratio", "equal_chance_feasible", "equal_chance_all_256", "P_best_abatement",
          "most_likely_has_best_abatement"]
with open(f"{OUT}/twolimit_qaoa_sweep.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f); wr.writerow(HEAD_S)
    for c in CASES:
        base = [c["elec"], f"{c['fB']:.2f}", f"{c['fS']:.2f}", z0(c["B"]), z0(c["S"]), c["nF"], plan_label(c["k_best"]),
                z0(c["abv"][c["k_best"]]), c["n_tied"], int(c["real_choice"]), int(c["dup_literal_of"] == c["id"]),
                c["same_F_as"], int(c["same_F_as"] == c["id"]), z0(c["P_fixed"]), int(c["fixed_sep"]), z0(c["P_need"]),
                z0(c["P_need_h"]), "not_run" if not c["real_choice"] else ("run" if c["P_min"] is not None else "not needed"),
                "" if math.isnan(c["min_viol"]) else z0(c["min_viol"])]
        if not c["rows"]:
            wr.writerow(base + ["not_run", "", "", "", "", "", "", "", z0(1.0 / c["nF"]), z0(1.0 / NPLAN), "", ""])
        for r in c["rows"]:
            wr.writerow(base + [r["method"], "" if r["penalty"] is None else z0(r["penalty"]), r["p"], z0(r["P_best"]),
                                z0(r["P_feasible"]), r["most_likely_plan"], r["most_likely_is_best"], z0(r["exp_ratio"]),
                                z0(1.0 / c["nF"]), z0(1.0 / NPLAN), z0(r["P_tied"]), r["ml_tied"]])

# summary. Two ways of counting: every real-choice case; and one case per (abatement table, feasible set).
# case_filter: "all cases run" = every case in which the method was run; "penalty separates" = only the cases in
# which, with the method's penalty, every infeasible plan has a higher energy than the best feasible plan
# (n_cases gives their number). O_min and O_hinge separate by construction wherever they are run, so their two
# rows are equal; for O_fixed the two rows differ. G has no penalty (one row).
HEAD_Y = ["counting", "method", "case_filter", "p", "n_cases", "mean_P_best", "min_P_best", "n_most_likely_is_best",
          "mean_P_feasible", "min_P_feasible", "mean_equal_chance_feasible", "mean_exp_abatement_ratio",
          "min_exp_abatement_ratio", "n_penalty_separates", "n_cases_with_tied_best", "mean_P_best_abatement",
          "min_P_best_abatement", "n_most_likely_has_best_abatement", "mean_equal_chance_best_abatement"]
SUMMARY = []
for counting, keep in [("all_real_choice_cases", lambda c: c["real_choice"]),
                       ("distinct_feasible_set", lambda c: c["real_choice"] and c["same_F_as"] == c["id"])]:
    for method in ("G", "O_fixed", "O_min", "O_hinge"):
        filters = [("all cases run", lambda r: True)]
        if method != "G":
            filters = [("penalty separates", lambda r: int(r["penalty_separates"]) == 1)] + filters
        for fname, fkeep in filters:
            for p in range(1, P_MAX_SWEEP + 1):
                sel = [(c, r) for c in CASES if keep(c) for r in c["rows"]
                       if r["method"] == method and r["p"] == p and fkeep(r)]
                if not sel:
                    continue
                pb = np.array([r["P_best"] for _, r in sel]); pf = np.array([r["P_feasible"] for _, r in sel])
                er = np.array([r["exp_ratio"] for _, r in sel])
                SUMMARY.append([counting, method, fname, p, len(sel), z0(pb.mean()), z0(pb.min()),
                                sum(r["most_likely_is_best"] for _, r in sel), z0(pf.mean()), z0(pf.min()),
                                z0(np.mean([1.0 / c["nF"] for c, _ in sel])), z0(er.mean()), z0(er.min()),
                                "" if method == "G" else sum(int(r["penalty_separates"]) for _, r in sel),
                                sum(c["n_tied"] > 1 for c, _ in sel), z0(np.mean([r["P_tied"] for _, r in sel])),
                                z0(min(r["P_tied"] for _, r in sel)), sum(r["ml_tied"] for _, r in sel),
                                z0(np.mean([c["n_tied"] / c["nF"] for c, _ in sel]))])
with open(f"{OUT}/twolimit_qaoa_summary.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f); wr.writerow(HEAD_Y); wr.writerows(SUMMARY)
print("\n  " + "  ".join(HEAD_Y))
for row in SUMMARY:
    print("  " + "  ".join(str(x) for x in row))
print("\n  BOTH COUNTS (quote together). 'single' = the single best plan (tie rule) is the most likely output;")
print("  'any best' = a plan with the best abatement is the most likely output (headline count).")
for row in SUMMARY:
    print(f"    {row[0]:22s} {row[1]:7s} {row[2]:17s} p={row[3]}: single {row[7]:3d} of {row[4]:3d} | any best "
          f"{row[17]:3d} of {row[4]:3d} | cases with tied best plans {row[14]}")

# =============================================================================
# PART 7 -- COMPARISON WITH THE EXACT CLASSICAL SOLVER (phase 7), IF ITS OUTPUT IS PRESENT
# =============================================================================
DEC = f"{OUT}/twolimit_decisions.csv"
if not os.path.exists(DEC):
    print(f"\nNOTE: {DEC} not found -> best plans of this script were NOT compared with the phase 7 solver.")
else:
    with open(DEC, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    cols = list(rows[0].keys()); low = {c.lower(): c for c in cols}
    def pick(*cands):
        for cand in cands:
            if cand in low:
                return low[cand]
        return None
    c_e = pick("electricity", "grid", "electricity_case"); c_b = pick("f_b", "fb", "f_budget", "budget_fraction")
    c_s = pick("f_s", "fs", "f_scrap", "scrap_fraction")
    c_p = pick("best_plan", "plan", "optimal_plan", "optimum_plan", "plan_optimal", "exact_plan", "brute_force_plan")
    if None in (c_e, c_b, c_s, c_p):
        print(f"\nWARNING: {DEC} exists but its columns {cols} were not recognised -> NOT compared.")
    else:
        mine = {(c["elec"], round(c["fB"], 2), round(c["fS"], 2)): plan_label(c["k_best"]) for c in CASES}
        n_cmp = 0
        for r in rows:
            key = (r[c_e], round(float(r[c_b]), 2), round(float(r[c_s]), 2))
            if key in mine:
                assert mine[key] == r[c_p].strip(), f"best plan differs from {DEC} at {key}: {mine[key]} vs {r[c_p]}"
                n_cmp += 1
        print(f"\nBest plans equal those in {DEC} in {n_cmp} of {len(mine)} cases compared "
              f"(columns {c_e}, {c_b}, {c_s}, {c_p}).")
        assert n_cmp == len(mine) == len(rows), "not every case of twolimit_decisions.csv was matched"
        if "n_tied_best" in cols:      # number of plans tied on the best abatement must also agree
            ties7 = {(r[c_e], round(float(r[c_b]), 2), round(float(r[c_s]), 2)): int(r["n_tied_best"]) for r in rows}
            for c in CASES:
                key = (c["elec"], round(c["fB"], 2), round(c["fS"], 2))
                assert ties7[key] == c["n_tied"], f"number of tied best plans differs from {DEC} at {key}"
            print(f"Number of tied best plans equal to {DEC} in all {len(CASES)} cases.")

with open(f"{OUT}/software_versions_phase8_twolimit_qaoa_colab.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["package", "version"])
    wr.writerows([["python", platform.python_version()], ["numpy", np.__version__], ["scipy", scipy.__version__]])

print(f"\nDone. Total run time {time.time() - T0:.1f} s. Outputs in {OUT}/: twolimit_qaoa_main.csv, "
      f"twolimit_qaoa_sweep.csv, twolimit_qaoa_summary.csv")
