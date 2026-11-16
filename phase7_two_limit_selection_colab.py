# -*- coding: utf-8 -*-
"""
phase7_two_limit_selection_colab.py
===================================
Plant choice with FOUR options per plant and TWO limits (capital in EUR, national scrap), solved exactly
by full enumeration, plus an uncertainty analysis on the Step III route model inputs (grid factor and
electricity need). No Step II value enters the plant choice.
No quantum code is in this file. Every result here is classical and exact.

HOW TO RUN IN GOOGLE COLAB
  Run phase2_route_estimates_colab.py, phase3_route_model_colab.py and
  phase6_mass_energy_balance_colab.py first (same runtime), or upload these files:
     outputs_phase3/plant_inputs.csv
     outputs_phase3/scenario_table.csv
     outputs_phase3/mass_energy_balance.csv
  Then run:   !python phase7_two_limit_selection_colab.py        (takes a few seconds)
Outputs (all in ./outputs_phase4/):
  twolimit_plant_options.csv        A, C, U for 4 plants x 4 options x 2 electricity cases
  twolimit_decisions.csv            exact best plan for 2 electricity cases x 10 capital limits x 5 scrap limits,
                                    with every plan tied with it on abatement (n_tied_best, tied_best_plans)
  twolimit_summary.csv              does the richer problem matter? counts per electricity case
  twolimit_uncertainty_plans.csv    each distinct best plan and its frequency over 2000 draws, 3 limit cases
  twolimit_uncertainty_summary.csv  share of draws in which the deterministic plan stays best; percentiles;
                                    one-input-fixed-at-a-time runs
  twolimit_flip_thresholds.csv      for 3 limit cases x 2 inputs (grid factor; electricity need at 0% scrap) x
                                    2 directions: the smallest change from the central value at which the best
                                    plan changes, the new plan, and the gap between best plan and runner-up
  software_versions_phase7_two_limit_selection_colab.csv   python and numpy versions of the run

DECISION
  Each plant i takes ONE option o in {N, S0, S50, S95}: do not convert, or convert to hydrogen direct
  reduction + electric arc furnace at scrap share s = 0, 0.50, 0.95.
  A plan is a 4-character string in plant order J, V, D, B with characters N, 0, 5, 9 (for example "N9N5").
    Abatement  A_i(o) = Q_i (GEI_i - I(s_o, g)) / 1e6        Mt CO2e/yr      ; A_i(N) = 0
    Capital    C_i(o) = Q_i capex(s_o) / 1e6                 million EUR     ; C_i(N) = 0
               capex(s) = 184 + (1 - s)(160 + 230)           EUR per t of yearly capacity
    Scrap      U_i(o) = Q_i m_scrap(s_o) / 1e6               Mt scrap/yr     ; U_i(N) = 0
               m_scrap(s) = scrap mass per t steel from the iron balance of phase 6 (stream
                            "scrap (iron-balance scaling M)" of outputs_phase3/mass_energy_balance.csv),
                            read for s = 0, 0.50 and 0.95; no straight line through zero is used
    Limits     sum C_i <= B = f_B x C_all0 ,  C_all0 = sum_i Q_i x 574 / 1e6
               sum U_i <= S = f_S x S_nat  ,  S_nat  = 34.2 Mt (India, 2024)
    Objective  maximise sum_i A_i(o_i); ties -> lower capital, then lower scrap.
    Ties       ALL feasible plans whose abatement equals the best abatement within 1e-9 Mt/yr are listed in
               tied_best_plans (separator "|"); best_plan is the one of them picked by the tie rule.
               With renewable electricity the route intensity does not depend on the scrap share, so many
               plans tie there.
    NO EPD VALUE IS USED IN THIS SCRIPT. No quantity multiplies a plant's abatement by EPD_i / GEI_i.
               Reason: an EPD (cradle to gate, product level, CO2e) and a GEI (plant level, per
               tonne of equivalent product) come from different accounting systems, and their ratio is below 1
               for two plants, so scaling one by the other gives values below the route intensity itself.
               This script does not read outputs_phase2/plant_epd_values.csv.
  Route model (the same as phase3_route_model_colab.py; constants restated below, nothing is imported):
    sec(s) = 3.48 - s (3.48 - 0.667) ;  I(s, g) = 0.053 + sec(s) g + 2.2 x 0.0561 + 0.140 g

RULES STATED HERE SO THAT THE PAPER CAN QUOTE THEM
  Exact solver: all 4^4 = 256 plans are listed; the feasible plan with the highest abatement is taken.
  Second method (check): the 256 plans are listed again as the integers 0..255, decoded as 8 bits
    (2 bits per plant, 00 = N, 01 = S0, 10 = S50, 11 = S95, plant J first) with a running comparison.
    The two methods must agree in every limit case (assert).
  Binding limit: a limit is called binding when the best abatement becomes strictly higher if that limit
    alone is removed and the other limit is kept. Label: "capital", "scrap", "both" (each one alone),
    "joint" (abatement rises only when both are removed together), "none".
    This rule is used in place of "the next step would break the limit" because, with four options per
    plant, a plan can leave unused capital or scrap that no single step can use.
  Greedy rule: start with no plant converted. A step changes the option of ONE plant. Among all steps
    that raise the total abatement and keep both limits, take the step with the highest extra abatement
    per extra EUR. A step that raises abatement and needs no extra capital (extra EUR <= 0) is taken
    before any step that needs capital (highest extra abatement first). Ties: higher extra abatement, then
    plant order J, V, D, B, then option order. Stop when no such step is left.
  Mixed plan: the converted plants of the best plan do not all use the same scrap share.
  GEI ranking: the best plan "follows GEI" when its k converted plants are the k plants with the highest GEI.

UNCERTAINTY (PART D)
  2000 draws, numpy default_rng(12345). Drawn independently: grid factor (uniform 0.710 to 0.727; the lower
  end is a reported value whose source document was not read), electricity need at 0% scrap (uniform within
  +/-5% of 3.48; no source), and the GEI of each plant (uniform between its 2025 and 2026 draft values).
  "Central value" in the one-input-fixed runs = the value of the deterministic case (0.727; 3.48; the 2026
  draft GEI). Note that 0.727 is the upper end of the grid-factor interval, not its middle.
  The fast (array) solver used for the draws is checked against the plain exact solver in ALL 2000 draws of
  every limit case (assert).
  WHAT PART D DOES NOT SHOW. The capital cost, the scrap mass per tonne of steel and both limits are NOT
  varied. The three ranges that are varied are narrow, and partly unsourced (the +/-5% has no source; 0.710 was
  not read in a primary document; the two gazette drafts are two points, not a range of error). A statement
  such as "the plan does not change in x% of the draws" is therefore NOT evidence that the plan is robust. It
  only says that these three small changes are too small to move it.

FLIP THRESHOLDS (PART E)
  For each of the three limit cases and for each of two inputs (grid factor; electricity need at 0% scrap),
  the input is moved away from its central value on a fine grid, upward and downward, with every other input
  at its central value, until the best plan changes. Reported: the first grid value at which the plan differs,
  the change in absolute and relative terms, the new plan, and the abatement gap between the best plan and the
  runner-up at the central values. Grid steps and search limits are study choices (constants FLIP_*). If the
  plan does not change inside the search limits, the row says so. These thresholds are computed, not assumed.

Tags: [VERIFIED] value read in the named source text (page marker of the extracted PDF text and a quote);
[ASSUMPTION] study choice to state in the paper; [UNVERIFIED] reported value, source document not read.
"""
import csv, os, itertools, platform, time
import numpy as np

T0 = time.time()

SEED = 12345
N_DRAWS = 2000
TOL = 1e-9

OUT = "outputs_phase4"
os.makedirs(OUT, exist_ok=True)
PLANTS_CSV = "outputs_phase3/plant_inputs.csv"
SCEN_CSV = "outputs_phase3/scenario_table.csv"
BAL_CSV = "outputs_phase3/mass_energy_balance.csv"
for need in (PLANTS_CSV, SCEN_CSV, BAL_CSV):
    if not os.path.exists(need):
        raise FileNotFoundError("Run phase 3 (route model) and phase 6 (mass and energy balance) first, "
                                "or upload the file: " + need)

# =============================================================================
# PART 1 -- CONSTANTS
# =============================================================================
# Vogl, Ahman & Nilsson (2018), J. Clean. Prod. 203, 736-745 (page = page of the PDF file)
SEC_H2DRI_0SCRAP = 3.48   # MWh/t liquid steel. PDF p. 4: "A specific energy consumption (SEC) of 3.48 MWh/t of
                          # liquid steel was calculated for the base case."                          [VERIFIED]
SEC_EAF_100SCRAP = 0.667  # MWh/t. PDF p. 5: "less energy (0.667 MWh/tLS) than from pure DRI"        [VERIFIED]
DIRECT_CO2 = 0.053        # t CO2/t. PDF p. 6: "graphite electrodes would result in emissions of
                          # 53 kg CO2 per" tonne of steel                                             [VERIFIED]
REST_FUEL_GJ = 2.2        # GJ/t. PDF p. 4: "For pelletizing and rolling an additional 2.2 GJ/t of fuel" [VERIFIED]
REST_ELEC_MWH = 0.140     # MWh/t. PDF p. 5: "140 kWh/t electricity are required irrespective of"     [VERIFIED]
# IPCC 2006 Guidelines, Vol. 2, Ch. 1, Table 1.4 (PDF p. 24): "Natural Gas 15.3 1 56 100 54 300 58 300"
EF_NATGAS = 0.0561        # t CO2/GJ                                                                  [VERIFIED]
                          # The 2.2 GJ fuel is taken as natural gas; Vogl does not name the fuel.     [ASSUMPTION]
# CEA, CO2 Baseline Database for the Indian Power Sector, User Guide Version 20.0 (PDF p. 4, Table S, and
# p. 15): "Weighted Average factor (Grid Emission Factor) = 1204.51 MT / 1655.7 BU = 0.727 tCO2/MWh (Table S)"
GRID_INDIA = 0.727        # t CO2/MWh, FY 2023-24                                                     [VERIFIED]
GRID_INDIA_V21_REPORTED = 0.710   # FY 2024-25; not checked in a CEA document when this script was run;
                                  # later confirmed by CEA User Guide Version 22.0 (Appendix C, Table B).
                                  # Used only as the lower end of the grid-factor interval in PART D.
GRID_RENEWABLE = 0.0      # t CO2/MWh, operation only                                                 [ASSUMPTION]
ELECTRICITY = [("India grid", GRID_INDIA), ("Renewable", GRID_RENEWABLE)]

# Capital cost. Vogl et al. (2018), PDF p. 6: "The CAPEX for a H-DR plant was calculated as 574 EUR per tonne
# capacity ... It includes 160 EUR/t capacity for the electrolyser, 230 EUR/t for the shaft and 184 EUR/t for
# the EAF". Same page: "Our CAPEX assessment is based on 2030 technology for the main components of the system."
# The paper gives no price year for the EUR values. European cost data, not Indian.
# The unit is "per tonne capacity" AS PRINTED in Vogl et al.; only the four numbers below are [VERIFIED].
CAPEX_ELECTROLYSER = 160.0   # EUR per tonne capacity                                                 [VERIFIED]
CAPEX_SHAFT = 230.0          # EUR per tonne capacity                                                 [VERIFIED]
CAPEX_EAF = 184.0            # EUR per tonne capacity                                                 [VERIFIED]
CAPEX_TOTAL_VOGL = 574.0     # EUR per tonne capacity                                                 [VERIFIED]
if CAPEX_ELECTROLYSER + CAPEX_SHAFT + CAPEX_EAF != CAPEX_TOTAL_VOGL:      # the three printed parts must add up
    raise ValueError("the three capital cost parts do not add up to the printed total of Vogl et al.")
# Reading "per tonne capacity" as per tonne of YEARLY STEEL capacity is our reading; the quoted sentence
# does not say "per year".                                                                          [ASSUMPTION]
# Using the plant output of 2023-24, which the gazette gives in tonnes of EQUIVALENT PRODUCT, as the steel
# capacity to be built is a second reading (equivalent product is not liquid steel).                [ASSUMPTION]
# Electrolyser and shaft capacity scale with the ore share (1 - s); the EAF does not.               [ASSUMPTION]
# No conversion to rupees, no cost of plant removal, no change of cost with plant size.             [ASSUMPTION]

# Scrap mass. Vogl et al. (2018), PDF p. 4: "If the EAF is charged with equal shares of HBI and scrap,
# 738 kg of pellets, 536 kg of scrap and 25 kg of hydrogen are needed per tonne of steel."
SCRAP_KG_PER_T_AT_50 = 536.0                             #                                           [VERIFIED]
# The 536 kg is used ONLY as a check. The scrap mass per tonne of steel is READ from the iron balance of
# phase6_mass_energy_balance_colab.py (PART 2 below). That balance uses the compositions printed by Vogl
# et al. (PDF p. 2: "Scrap charged to the EAF contains 95% iron and 5% inert substances."; PDF p. 3: "The
# liquid steel product contains only iron") and closes the iron balance at 1000 kg iron per tonne of steel
# at every scrap share (by construction of that scaling). Phase 6 also lists a second scaling ("proportional scaling P",
# 1072 kg x s); it does not close the iron balance away from s = 0.50 and is NOT used here.
# In phase 6 the scrap share s is read as scrap mass / (scrap mass + HBI mass) of the furnace charge. [ASSUMPTION]
SCRAP_STREAM = "scrap (iron-balance scaling M)"          # stream name as written by phase 6
SCRAP_IRON_STREAM = "iron in scrap (M)"
X_FE_SCRAP = 0.95                                        # iron share of scrap, Vogl PDF p. 2 (quote above) [VERIFIED]

# National scrap. EY-Parthenon, WWF-India & CII-Green Business Centre (February 2026), "Closing the loop",
# PDF p. 18: "In 2024, India's total ferrous scrap consumption stood at 34.2 million tons." The sentence before
# it reads: "India's ferrous scrap consumption in steelmaking rose by 10% year-on-year in 2024, reaching above
# 34 MMT." The four regional values on the same page add up to 34.2 (9.4 + 9.2 + 7.8 + 7.8). Source named
# there: Bigmint, January 2025. Fig. 1 (PDF p. 10) gives steelmaking scrap as BOF 6 + EAF 7.5 + IF 21 MMT
# (= 34.5) and foundries 6.50 MMT, total 41 MMT. So the 34.2 Mt covers steelmaking in all three furnace types
# (induction furnace, electric arc furnace, basic oxygen furnace) and NOT foundries. The report does not say
# this in one sentence; it is read from the two places together.
S_NAT_MT = 34.2           # Mt ferrous scrap consumed, India, 2024                                    [VERIFIED]
# The scrap limit S = f_S x S_NAT_MT is a share of scrap that other steelmakers use TODAY. It is a size
# reference, not a forecast of the scrap that the four plants could buy.                             [ASSUMPTION]

OPTIONS = [("N", None), ("0", 0.0), ("5", 0.50), ("9", 0.95)]     # code, scrap share        [ASSUMPTION: 3 shares]
OPT_NAME = {"N": "not converted", "0": "S0", "5": "S50", "9": "S95"}
F_B_GRID = [round(0.10 * k, 2) for k in range(1, 11)]
F_S_GRID = [0.10, 0.25, 0.50, 0.75, 1.00]
MAIN = ("India grid", 0.50, 0.25)                                 # main case                 [ASSUMPTION]
UNC_CASES = [(0.50, 0.25), (0.50, 0.50), (0.50, 1.00)]            # limit cases of PART D     [ASSUMPTION]
SEC_REL_HALF_WIDTH = 0.05   # electricity need at 0% scrap: uniform within +/-5% of 3.48. NO SOURCE.  [ASSUMPTION]
# PART E, flip thresholds: grid step and search limits of the one-input scans.                      [ASSUMPTION]
FLIP_G_STEP = 0.0001        # t CO2/MWh; grid factor scanned from 0 to 2 x 0.727
FLIP_G_MAX_FACTOR = 2.0
FLIP_SEC_STEP = 0.001       # MWh/t; electricity need at 0% scrap scanned from 0.667 (the pure-scrap value, below
FLIP_SEC_MAX_FACTOR = 2.0   # which the straight line has no meaning) to 2 x 3.48

# =============================================================================
# PART 2 -- INPUTS AND MODEL
# =============================================================================
with open(PLANTS_CSV, encoding="utf-8") as f:
    _rows = list(csv.DictReader(f))
PLANT_NAMES = [r["plant"] for r in _rows]
for name, key in zip(PLANT_NAMES, ["Jamshedpur", "Vijayanagar", "Dolvi", "Bhilai"]):
    assert key in name, "plant order in plant_inputs.csv is not J, V, D, B"
SHORT = ["J", "V", "D", "B"]
Q = [float(r["output_2023_24_t_rev2026"]) for r in _rows]          # t/yr, 2026 revised draft
GEI_2026 = [float(r["GEI_2023_24_rev2026"]) for r in _rows]        # t CO2e/t, 2026 revised draft (main)
GEI_2025 = [float(r["GEI_2023_24_draft2025"]) for r in _rows]      # t CO2e/t, 2025 draft
# The four 2025-draft values 2.3804, 2.6293, 2.6662, 3.1487 were also seen in MoEFCC_2025_405E (PDF p. 72)
# and the 2026 rows (output, GEI) in MoEFCC_2026_517E, next to the registration numbers of the CSV.
NP = len(Q)
assert NP == 4

for _n, _q, _g6, _g5 in zip(PLANT_NAMES, Q, GEI_2026, GEI_2025):
    if not (_q > 0 and _g6 > 0 and _g5 > 0):
        raise ValueError("plant without a positive output or GEI in " + PLANTS_CSV + ": " + _n)


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


SCRAP_T_PER_T = read_scrap_t_per_t(BAL_CSV)
# Note: at s = 0.95 the scrap MASS (about 1002 kg) is above 1000 kg of steel because scrap is 95 % iron; the
# IRON in that scrap (about 952 kg) is below the 1000 kg of iron in the steel, which is what is asserted.


def sec(s, sec0=SEC_H2DRI_0SCRAP):                   # MWh/t; linear between the two end points
    return sec0 - s * (sec0 - SEC_EAF_100SCRAP)


def intensity(s, g, sec0=SEC_H2DRI_0SCRAP):          # t CO2 / t liquid steel
    return DIRECT_CO2 + sec(s, sec0) * g + REST_FUEL_GJ * EF_NATGAS + REST_ELEC_MWH * g


def capex(s):                                        # EUR per t of yearly capacity
    return CAPEX_EAF + (1.0 - s) * (CAPEX_ELECTROLYSER + CAPEX_SHAFT)


def m_scrap(s):                                      # t scrap per t steel (iron balance of phase 6)
    return SCRAP_T_PER_T[round(s, 2)]


# Link to Phase III: the restated constants must give the Phase III table values (6 decimals).
with open(SCEN_CSV, encoding="utf-8") as f:
    _sc = {(float(r["scrap_share"]), r["electricity"]): float(r["I_new_tCO2_per_t"]) for r in csv.DictReader(f)}
for _s in (0.95, 0.50, 0.0):
    for _e, _g in ELECTRICITY:
        assert round(intensity(_s, _g), 6) == round(_sc[(_s, _e)], 6), \
            f"I({_s}, {_e}) differs from outputs_phase3/scenario_table.csv"
assert abs(capex(0.0) - CAPEX_TOTAL_VOGL) < 1e-12

C_ALL0 = sum(q * CAPEX_TOTAL_VOGL for q in Q) / 1e6      # million EUR: all four plants at 0% scrap


def tables(g, gei, sec0=SEC_H2DRI_0SCRAP):
    """A, C, U as [plant][option] lists. Option order N, S0, S50, S95."""
    A, C, U = [], [], []
    for i in range(NP):
        a, c, u = [0.0], [0.0], [0.0]
        for _, s in OPTIONS[1:]:
            a.append(Q[i] * (gei[i] - intensity(s, g, sec0)) / 1e6)
            c.append(Q[i] * capex(s) / 1e6)
            u.append(Q[i] * m_scrap(s) / 1e6)
        A.append(a); C.append(c); U.append(u)
    return A, C, U


def plan_str(p):
    return "".join(OPTIONS[o][0] for o in p)


def plan_words(p):
    parts = [f"{SHORT[i]}:{OPT_NAME[OPTIONS[o][0]]}" for i, o in enumerate(p) if o != 0]
    return " + ".join(parts) if parts else "none"


def totals(p, A, C, U):
    return (sum(A[i][p[i]] for i in range(NP)), sum(C[i][p[i]] for i in range(NP)),
            sum(U[i][p[i]] for i in range(NP)))


ALL_PLANS = list(itertools.product(range(4), repeat=NP))     # 256 plans, J first


def solve(A, C, U, B=None, S=None):
    """Exact: best feasible plan. B or S = None means that limit is ignored. Returns plan, a, c, u, n_feasible."""
    best, key_best, n = None, None, 0
    for p in ALL_PLANS:
        a, c, u = totals(p, A, C, U)
        if (B is not None and c > B + TOL) or (S is not None and u > S + TOL):
            continue
        n += 1
        key = (-round(a, 9), round(c, 9), round(u, 9))
        if key_best is None or key < key_best:
            best, key_best = p, key
    return (best,) + totals(best, A, C, U) + (n,)


def tied_best(A, C, U, B, S, a_best):
    """All feasible plans whose abatement equals a_best within TOL, in the order of the tie rule."""
    out = []
    for p in ALL_PLANS:
        a, c, u = totals(p, A, C, U)
        if c > B + TOL or u > S + TOL or abs(a - a_best) > TOL:
            continue
        out.append((round(c, 9), round(u, 9), p))
    return [p for _, _, p in sorted(out)]


def solve_bits(A, C, U, B, S):
    """Second method, written differently: integers 255..0 decoded as 8 bits, running comparison."""
    best = None
    for k in range(255, -1, -1):
        bits = format(k, "08b")
        p = tuple(int(bits[2 * i:2 * i + 2], 2) for i in range(NP))
        a = c = u = 0.0
        for i in range(NP - 1, -1, -1):
            a += A[i][p[i]]; c += C[i][p[i]]; u += U[i][p[i]]
        if c > B + TOL or u > S + TOL:
            continue
        if best is None:
            best = (p, a, c, u); continue
        _, a0, c0, u0 = best
        if a > a0 + TOL:
            best = (p, a, c, u)
        elif abs(a - a0) <= TOL:
            if c < c0 - TOL or (abs(c - c0) <= TOL and u < u0 - TOL):
                best = (p, a, c, u)
    return best[0]


def greedy(A, C, U, B, S):
    """Greedy rule of the docstring."""
    p = [0] * NP
    while True:
        a0, c0, u0 = totals(p, A, C, U)
        cand = []
        for i in range(NP):
            for o in range(4):
                if o == p[i]:
                    continue
                dA, dC, dU = A[i][o] - A[i][p[i]], C[i][o] - C[i][p[i]], U[i][o] - U[i][p[i]]
                if dA <= TOL or c0 + dC > B + TOL or u0 + dU > S + TOL:
                    continue
                free = dC <= TOL
                ratio = float("inf") if free else dA / dC
                cand.append((0 if free else 1, -ratio, -dA, i, o))
        if not cand:
            return tuple(p)
        _, _, _, i, o = min(cand)
        p[i] = o


def binding_label(a_best, a_noB, a_noS, a_free):
    cap, scr = a_noB > a_best + TOL, a_noS > a_best + TOL
    if cap and scr:
        return "both"
    if cap:
        return "capital"
    if scr:
        return "scrap"
    return "joint" if a_free > a_best + TOL else "none"


def follows_gei(p, gei):
    conv = {i for i in range(NP) if p[i] != 0}
    order = sorted(range(NP), key=lambda i: -gei[i])
    return conv == set(order[:len(conv)])


def is_mixed(p):
    return len({o for o in p if o != 0}) > 1


# =============================================================================
# PART A -- OPTION TABLES
# =============================================================================
print(f"C_all0 = {C_ALL0:.3f} million EUR (all four plants at 0% scrap); S_nat = {S_NAT_MT} Mt")
print(f"capex(s): 0 -> {capex(0):.1f}, 0.50 -> {capex(0.5):.1f}, 0.95 -> {capex(0.95):.1f} EUR/t")
print(f"scrap t per t steel (iron balance, phase 6): s = 0 -> {m_scrap(0):.6f}, 0.50 -> {m_scrap(0.5):.6f}, "
      f"0.95 -> {m_scrap(0.95):.6f}  (Vogl prints 0.536 at s = 0.50: "
      f"{100 * (m_scrap(0.5) / 0.536 - 1):+.3f} %)")
with open(f"{OUT}/twolimit_plant_options.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["electricity", "grid_tCO2_per_MWh", "plant", "plant_code", "option", "scrap_share",
                 "output_t_per_yr", "GEI_tCO2e_per_t", "I_new_tCO2_per_t", "capex_EUR_per_t", "scrap_t_per_t",
                 "abatement_Mt_per_yr", "capital_million_EUR", "scrap_Mt_per_yr",
                 "abatement_t_per_thousand_EUR", "abatement_t_per_t_scrap"])
    for ename, g in ELECTRICITY:
        A, C, U = tables(g, GEI_2026)
        print(f"\nOption table, {ename} (A Mt CO2e/yr | C million EUR | U Mt scrap/yr)")
        for i in range(NP):
            line = f"  {SHORT[i]}"
            for o, (code, s) in enumerate(OPTIONS):
                if o == 0:
                    wr.writerow([ename, f"{g:.6f}", PLANT_NAMES[i], SHORT[i], "N", "", int(Q[i]),
                                 f"{GEI_2026[i]:.6f}", "", f"{0:.6f}", f"{0:.6f}", f"{0:.6f}", f"{0:.6f}",
                                 f"{0:.6f}", "", ""])
                    continue
                per_eur = A[i][o] * 1e6 / (C[i][o] * 1e3)                       # t CO2e per thousand EUR
                per_scrap = "" if U[i][o] == 0 else f"{A[i][o] / U[i][o]:.6f}"
                wr.writerow([ename, f"{g:.6f}", PLANT_NAMES[i], SHORT[i], "S" + str(int(round(s * 100))), s,
                             int(Q[i]), f"{GEI_2026[i]:.6f}", f"{intensity(s, g):.6f}", f"{capex(s):.6f}",
                             f"{m_scrap(s):.6f}", f"{A[i][o]:.6f}", f"{C[i][o]:.6f}", f"{U[i][o]:.6f}",
                             f"{per_eur:.6f}", per_scrap])
                line += f" | {OPT_NAME[code]}: {A[i][o]:7.3f} {C[i][o]:8.1f} {U[i][o]:6.3f}"
            print(line)

# =============================================================================
# PART B -- EXACT SOLUTION FOR EVERY LIMIT CASE
# =============================================================================
DEC = []
for ename, g in ELECTRICITY:
    A, C, U = tables(g, GEI_2026)
    p_free, a_free, c_free, u_free, _ = solve(A, C, U)
    for fB in F_B_GRID:
        for fS in F_S_GRID:
            B, S = fB * C_ALL0, fS * S_NAT_MT
            p, a, c, u, nfeas = solve(A, C, U, B, S)
            assert solve_bits(A, C, U, B, S) == p, f"second method differs: {ename} f_B={fB} f_S={fS}"
            assert c <= B + TOL and u <= S + TOL
            ties = tied_best(A, C, U, B, S, a)
            assert ties[0] == p and 1 <= len(ties) <= nfeas, "tie list does not start with the best plan"
            p_noS, a_noS, c_noS, u_noS, _ = solve(A, C, U, B, None)
            p_noB, a_noB, c_noB, u_noB, _ = solve(A, C, U, None, S)
            assert a_noS >= a - TOL and a_noB >= a - TOL and a_free >= max(a_noS, a_noB) - TOL
            pg = greedy(A, C, U, B, S)
            ag, cg, ug = totals(pg, A, C, U)
            assert cg <= B + TOL and ug <= S + TOL and ag <= a + TOL
            short = 0.0 if a <= TOL else (a - ag) / a * 100
            DEC.append(dict(e=ename, g=g, fB=fB, fS=fS, B=B, S=S, p=p, a=a, c=c, u=u, nfeas=nfeas,
                            bind=binding_label(a, a_noB, a_noS, a_free), pg=pg, ag=ag, cg=cg, ug=ug, short=short,
                            p_noS=p_noS, a_noS=a_noS, c_noS=c_noS, u_noS=u_noS,
                            p_noB=p_noB, a_noB=a_noB, c_noB=c_noB, u_noB=u_noB,
                            mixed=is_mixed(p), gei_ok=follows_gei(p, GEI_2026), main=(ename, fB, fS) == MAIN,
                            ties=ties))

with open(f"{OUT}/twolimit_decisions.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["electricity", "f_B", "f_S", "capital_limit_million_EUR", "scrap_limit_Mt", "main_case",
                 "n_feasible_plans_of_256", "best_plan", "best_plan_words", "abatement_Mt_per_yr",
                 "capital_used_million_EUR", "scrap_used_Mt", "capital_used_share_of_limit",
                 "scrap_used_share_of_limit", "binding_limit", "best_plan_mixes_scrap_shares",
                 "converted_plants_are_top_GEI", "greedy_plan", "greedy_abatement_Mt_per_yr",
                 "greedy_shortfall_pct", "greedy_is_optimal",
                 "best_plan_no_scrap_limit", "abatement_no_scrap_limit", "scrap_needed_no_scrap_limit_Mt",
                 "scrap_limit_changes_plan", "abatement_lost_to_scrap_limit_Mt",
                 "best_plan_no_capital_limit", "abatement_no_capital_limit", "capital_needed_no_capital_limit",
                 "capital_limit_changes_plan", "abatement_lost_to_capital_limit_Mt",
                 "n_tied_best", "tied_best_plans"])
    for d in DEC:
        wr.writerow([d["e"], f"{d['fB']:.2f}", f"{d['fS']:.2f}", f"{d['B']:.6f}", f"{d['S']:.6f}", int(d["main"]),
                     d["nfeas"], plan_str(d["p"]), plan_words(d["p"]), f"{d['a']:.6f}", f"{d['c']:.6f}",
                     f"{d['u']:.6f}", f"{d['c'] / d['B']:.6f}", f"{d['u'] / d['S']:.6f}", d["bind"],
                     int(d["mixed"]), int(d["gei_ok"]), plan_str(d["pg"]), f"{d['ag']:.6f}",
                     f"{d['short']:.6f}", int(d["pg"] == d["p"] or d["short"] <= 1e-9),
                     plan_str(d["p_noS"]), f"{d['a_noS']:.6f}", f"{d['u_noS']:.6f}", int(d["p_noS"] != d["p"]),
                     f"{d['a_noS'] - d['a']:.6f}",
                     plan_str(d["p_noB"]), f"{d['a_noB']:.6f}", f"{d['c_noB']:.6f}", int(d["p_noB"] != d["p"]),
                     f"{d['a_noB'] - d['a']:.6f}",
                     len(d["ties"]), "|".join(plan_str(q) for q in d["ties"])])


def show(rows, title):
    print("\n" + title)
    print("  f_B  f_S | best  abat.   capital  scrap | binding | feas | greedy short% | noScrapLim  noCapLim"
          " | tied")
    for d in rows:
        print(f"  {d['fB']:.2f} {d['fS']:.2f} | {plan_str(d['p'])} {d['a']:7.3f} {d['c']:8.1f} {d['u']:6.3f} | "
              f"{d['bind']:7s} | {d['nfeas']:4d} | {plan_str(d['pg'])} {d['short']:6.2f} | "
              f"{plan_str(d['p_noS'])} {d['a_noS']:6.3f}  {plan_str(d['p_noB'])} {d['a_noB']:6.3f} | "
              f"{len(d['ties']):4d}")


show([d for d in DEC if d["main"]], f"MAIN CASE {MAIN}")
show([d for d in DEC if d["e"] == "India grid"], "ALL INDIA-GRID CASES (plan order J V D B; N, 0, 5, 9)")
show([d for d in DEC if d["e"] == "Renewable"], "ALL RENEWABLE CASES")

# =============================================================================
# PART C -- DOES THE RICHER PROBLEM MATTER?
# =============================================================================
with open(f"{OUT}/twolimit_summary.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["electricity", "n_limit_cases", "n_distinct_best_plans", "n_best_plan_converts_nothing",
                 "n_best_plan_mixes_scrap_shares", "n_best_plan_uses_S0", "n_best_plan_uses_S50",
                 "n_best_plan_uses_S95", "n_scrap_limit_changes_plan", "n_scrap_limit_lowers_abatement",
                 "n_capital_limit_changes_plan", "n_capital_limit_lowers_abatement",
                 "n_binding_capital_only", "n_binding_scrap_only", "n_binding_both", "n_binding_joint",
                 "n_binding_none", "n_greedy_not_optimal", "greedy_largest_shortfall_pct",
                 "greedy_mean_shortfall_pct_all_cases", "greedy_mean_shortfall_pct_nonoptimal_cases",
                 "n_converted_plants_not_top_GEI", "n_cases_with_tied_best_plans", "largest_n_tied_best"])
    print("\nPART C -- counts over the 50 limit cases per electricity case")
    for ename, _ in ELECTRICITY:
        R = [d for d in DEC if d["e"] == ename]
        bad = [d["short"] for d in R if d["short"] > 1e-9]
        nb = lambda lab: sum(d["bind"] == lab for d in R)
        row = [ename, len(R), len({d["p"] for d in R}), sum(all(o == 0 for o in d["p"]) for d in R),
               sum(d["mixed"] for d in R), sum(1 in d["p"] for d in R), sum(2 in d["p"] for d in R),
               sum(3 in d["p"] for d in R), sum(d["p_noS"] != d["p"] for d in R),
               sum(d["a_noS"] > d["a"] + TOL for d in R), sum(d["p_noB"] != d["p"] for d in R),
               sum(d["a_noB"] > d["a"] + TOL for d in R), nb("capital"), nb("scrap"), nb("both"), nb("joint"),
               nb("none"), len(bad), f"{max(bad) if bad else 0:.6f}", f"{sum(bad) / len(R):.6f}",
               f"{(sum(bad) / len(bad)) if bad else 0:.6f}", sum(not d["gei_ok"] for d in R),
               sum(len(d["ties"]) > 1 for d in R), max(len(d["ties"]) for d in R)]
        wr.writerow(row)
        print(f"  {ename}: distinct best plans {row[2]}; mixed scrap shares {row[4]}/50; uses S0 {row[5]}, "
              f"S50 {row[6]}, S95 {row[7]}; scrap limit changes plan {row[8]} (lowers abatement {row[9]}); "
              f"capital limit changes plan {row[10]} (lowers abatement {row[11]}); binding capital/scrap/both/"
              f"joint/none {row[12]}/{row[13]}/{row[14]}/{row[15]}/{row[16]}; greedy not optimal {row[17]}, "
              f"largest {float(row[18]):.2f}%, mean {float(row[19]):.2f}% (all), {float(row[20]):.2f}% "
              f"(non-optimal); converted plants not top-GEI {row[21]}; cases with tied best plans {row[22]} "
              f"(largest number of tied plans {row[23]})")

# =============================================================================
# PART D -- UNCERTAINTY (India grid; three limit cases)
# =============================================================================
rng = np.random.default_rng(SEED)
D_G = rng.uniform(GRID_INDIA_V21_REPORTED, GRID_INDIA, N_DRAWS)
D_SEC = rng.uniform(SEC_H2DRI_0SCRAP * (1 - SEC_REL_HALF_WIDTH), SEC_H2DRI_0SCRAP * (1 + SEC_REL_HALF_WIDTH),
                    N_DRAWS)
_lo = np.minimum(GEI_2025, GEI_2026); _hi = np.maximum(GEI_2025, GEI_2026)
D_GEI = rng.uniform(_lo, _hi, (N_DRAWS, NP))          # a plant with the same value in both drafts stays fixed

Q_ARR = np.array(Q)
S_OPT = np.array([s for _, s in OPTIONS[1:]])
PLANS_ARR = np.array(ALL_PLANS)
_, C_det, U_det = tables(GRID_INDIA, GEI_2026)
C_TOT = np.array([sum(C_det[i][p[i]] for i in range(NP)) for p in ALL_PLANS])
U_TOT = np.array([sum(U_det[i][p[i]] for i in range(NP)) for p in ALL_PLANS])
ORDER = np.lexsort((np.round(U_TOT, 9), np.round(C_TOT, 9)))     # by capital, then scrap (tie rule)


def abatement_all_plans(g, sec0, gei):
    """Total abatement of all 256 plans for every draw: array (n_draws, 256)."""
    sec_s = sec0[:, None] - S_OPT[None, :] * (sec0[:, None] - SEC_EAF_100SCRAP)
    I = DIRECT_CO2 + sec_s * g[:, None] + REST_FUEL_GJ * EF_NATGAS + REST_ELEC_MWH * g[:, None]   # (n, 3)
    A = np.zeros((len(g), NP, 4))
    A[:, :, 1:] = Q_ARR[None, :, None] * (gei[:, :, None] - I[:, None, :]) / 1e6
    tot = np.zeros((len(g), len(ALL_PLANS)))
    for i in range(NP):
        tot += A[:, i, PLANS_ARR[:, i]]
    return tot


def best_per_draw(tot, B, S):
    """Index of the best feasible plan per draw (ties: lower capital, then lower scrap)."""
    feas = (C_TOT <= B + TOL) & (U_TOT <= S + TOL)
    idx = ORDER[feas[ORDER]]
    sub = tot[:, idx]
    first = np.argmax(sub >= sub.max(axis=1, keepdims=True) - TOL, axis=1)
    return idx[first]


FIX_SETS = [("none (all inputs drawn)", set()), ("grid factor", {"g"}), ("electricity need at 0% scrap", {"sec"}),
            ("GEI of all four plants", {"gei0", "gei1", "gei2", "gei3"}), ("GEI Jamshedpur", {"gei0"}),
            ("GEI Vijayanagar", {"gei1"}), ("GEI Dolvi", {"gei2"}), ("GEI Bhilai", {"gei3"})]


def inputs(fixed):
    g = np.full(N_DRAWS, GRID_INDIA) if "g" in fixed else D_G
    s0 = np.full(N_DRAWS, SEC_H2DRI_0SCRAP) if "sec" in fixed else D_SEC
    gei = D_GEI.copy()
    for i in range(NP):
        if f"gei{i}" in fixed:
            gei[:, i] = GEI_2026[i]
    return g, s0, gei


pct = lambda x: [float(v) for v in np.percentile(x, [5, 50, 95])]
plan_rows, sum_rows = [], []
print(f"\nPART D -- {N_DRAWS} draws, India grid. Inputs: grid {GRID_INDIA_V21_REPORTED}-{GRID_INDIA}; "
      f"SEC0 +/-{SEC_REL_HALF_WIDTH:.0%}; GEI between the 2025 and 2026 drafts")
print("  NOT VARIED: capital cost, scrap mass per tonne of steel, and both limits. The three ranges that are varied "
      "are narrow and partly unsourced.")
print("  Therefore 'the plan does not change' in these draws is NOT evidence of robustness. See PART E for the size "
      "of change that does move the plan.")
N_EXACT_CHECKED = 0
for fB, fS in UNC_CASES:
    B, S = fB * C_ALL0, fS * S_NAT_MT
    det = next(d for d in DEC if d["e"] == "India grid" and d["fB"] == fB and d["fS"] == fS)
    det_idx = ALL_PLANS.index(det["p"])
    case = f"India grid, f_B={fB:.2f}, f_S={fS:.2f}"
    res = {}
    for label, fixed in FIX_SETS:
        g, s0, gei = inputs(fixed)
        tot = abatement_all_plans(g, s0, gei)
        bi = best_per_draw(tot, B, S)
        a_best = tot[np.arange(N_DRAWS), bi]
        res[label] = dict(bi=bi, a=a_best, a_det=tot[:, det_idx], share=float(np.mean(bi == det_idx)))
        if not fixed:      # check of the fast solver against the plain exact solver in ALL draws
            for n in range(N_DRAWS):
                A1, C1, U1 = tables(float(g[n]), list(gei[n]), float(s0[n]))
                assert solve(A1, C1, U1, B, S)[0] == ALL_PLANS[bi[n]], f"fast solver differs from exact solver (draw {n})"
                N_EXACT_CHECKED += 1
            one = np.ones(1)
            t0 = abatement_all_plans(one * GRID_INDIA, one * SEC_H2DRI_0SCRAP, np.array([GEI_2026]))
            assert best_per_draw(t0, B, S)[0] == det_idx and abs(t0[0, det_idx] - det["a"]) < 1e-9
    base = res["none (all inputs drawn)"]
    # frequency of each distinct best plan
    ids, counts = np.unique(base["bi"], return_counts=True)
    for k in np.argsort(-counts):
        p = ALL_PLANS[ids[k]]
        plan_rows.append([case, f"{fB:.2f}", f"{fS:.2f}", plan_str(p), plan_words(p), int(counts[k]),
                          f"{counts[k] / N_DRAWS:.6f}", int(ids[k] == det_idx), f"{C_TOT[ids[k]]:.6f}",
                          f"{U_TOT[ids[k]]:.6f}", f"{float(np.mean(base['a'][base['bi'] == ids[k]])):.6f}"])
    # sensitivity: which single input, when fixed, changes the result most
    cand = [l for l, fx in FIX_SETS if fx and l != "GEI of all four plants"]
    by_share = max(cand, key=lambda l: abs(res[l]["share"] - base["share"]))
    by_share = by_share if abs(res[by_share]["share"] - base["share"]) > 0 else "no single input changes the share"
    by_sd = min(cand, key=lambda l: float(np.std(res[l]["a"])))
    print(f"\n  {case}: deterministic best plan {plan_str(det['p'])} ({det['a']:.3f} Mt/yr)")
    for label, fixed in FIX_SETS:
        r = res[label]
        regret = r["a"] - r["a_det"]
        pa, pd_ = pct(r["a"]), pct(r["a_det"])
        sum_rows.append([case, f"{fB:.2f}", f"{fS:.2f}", label, N_DRAWS, plan_str(det["p"]), f"{det['a']:.6f}",
                         f"{r['share']:.6f}", len(np.unique(r["bi"])),
                         f"{pa[0]:.6f}", f"{pa[1]:.6f}", f"{pa[2]:.6f}", f"{float(np.std(r['a'])):.6f}",
                         f"{pd_[0]:.6f}", f"{pd_[1]:.6f}", f"{pd_[2]:.6f}",
                         f"{float(regret.mean()):.6f}", f"{float(regret.max()):.6f}", by_share, by_sd])
        print(f"    fixed: {label:30s} share det. plan best {r['share']:.4f} | plans {len(np.unique(r['bi']))} | "
              f"abat. p5/p50/p95 {pa[0]:.3f}/{pa[1]:.3f}/{pa[2]:.3f} sd {np.std(r['a']):.3f} | "
              f"loss if det. plan kept: mean {regret.mean():.4f} max {regret.max():.4f}")
    print(f"    most sensitive input by share of draws: {by_share}; by spread of abatement: {by_sd}")

with open(f"{OUT}/twolimit_uncertainty_plans.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["limit_case", "f_B", "f_S", "best_plan", "best_plan_words", "n_draws_best", "share_of_draws",
                 "is_deterministic_best_plan", "capital_million_EUR", "scrap_Mt",
                 "mean_abatement_when_best_Mt_per_yr"])
    wr.writerows(plan_rows)
with open(f"{OUT}/twolimit_uncertainty_summary.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["limit_case", "f_B", "f_S", "input_fixed_at_deterministic_value", "n_draws",
                 "deterministic_best_plan", "deterministic_abatement_Mt_per_yr",
                 "share_draws_deterministic_plan_is_best", "n_distinct_best_plans",
                 "abatement_p5", "abatement_p50", "abatement_p95", "abatement_sd",
                 "deterministic_plan_abatement_p5", "deterministic_plan_abatement_p50",
                 "deterministic_plan_abatement_p95", "mean_abatement_lost_if_deterministic_plan_kept",
                 "max_abatement_lost_if_deterministic_plan_kept",
                 "most_sensitive_single_input_by_plan_share", "most_sensitive_single_input_by_abatement_sd"])
    wr.writerows(sum_rows)

print(f"\n  fast solver = exact solver in all {N_EXACT_CHECKED} draws checked "
      f"({len(UNC_CASES)} limit cases x {N_DRAWS} draws)")

# =============================================================================
# PART E -- FLIP THRESHOLDS: how far must ONE input move before the best plan changes?
#   One input is moved on a fine grid away from its central value; all other inputs stay central (grid 0.727,
#   electricity need 3.48, GEI of the 2026 draft). The scan uses the fast solver; the reported grid value and the
#   grid value before it are then solved again with the plain exact solver (assert).
# =============================================================================
def scan_values(central, step, lo, hi, direction):
    """grid values central + k x step (direction 'up') or central - k x step ('down'), k = 1, 2, ..., inside [lo, hi]"""
    n = int(np.floor(((hi - central) if direction == "up" else (central - lo)) / step + 1e-9))
    k = np.arange(1, n + 1)
    return central + k * step if direction == "up" else central - k * step


FLIP_INPUTS = [
    ("grid factor", "t CO2/MWh", GRID_INDIA, FLIP_G_STEP, 0.0, FLIP_G_MAX_FACTOR * GRID_INDIA),
    ("electricity need at 0% scrap", "MWh/t", SEC_H2DRI_0SCRAP, FLIP_SEC_STEP, SEC_EAF_100SCRAP,
     FLIP_SEC_MAX_FACTOR * SEC_H2DRI_0SCRAP),
]
flip_rows = []
print("\nPART E -- flip thresholds (one input moved, the others at central values; India grid case)")
for fB, fS in UNC_CASES:
    B, S = fB * C_ALL0, fS * S_NAT_MT
    A0, C0, U0 = tables(GRID_INDIA, GEI_2026)
    p0 = solve(A0, C0, U0, B, S)[0]
    # runner-up at the central values: second feasible plan in the order of the solver (abatement, capital, scrap)
    ranked = sorted((-round(totals(p, A0, C0, U0)[0], 9), round(totals(p, A0, C0, U0)[1], 9),
                     round(totals(p, A0, C0, U0)[2], 9), p) for p in ALL_PLANS
                    if totals(p, A0, C0, U0)[1] <= B + TOL and totals(p, A0, C0, U0)[2] <= S + TOL)
    assert ranked[0][3] == p0, "ranking and solver disagree on the best plan"
    a0 = totals(p0, A0, C0, U0)[0]
    if len(ranked) > 1:
        p_ru = ranked[1][3]; a_ru = totals(p_ru, A0, C0, U0)[0]
        ru_txt, a_ru_txt, gap_txt, gap_pct_txt = plan_str(p_ru), f"{a_ru:.6f}", f"{a0 - a_ru:.6f}", f"{100 * (a0 - a_ru) / a0:.6f}"
    else:
        ru_txt, a_ru_txt, gap_txt, gap_pct_txt = "none (one feasible plan)", "", "", ""
    case = f"India grid, f_B={fB:.2f}, f_S={fS:.2f}"
    print(f"  {case}: best {plan_str(p0)} {a0:.6f} Mt/yr; runner-up {ru_txt} {a_ru_txt}; gap {gap_txt} Mt/yr ({gap_pct_txt} %)")
    for name, unit, central, step, lo, hi in FLIP_INPUTS:
        for direction in ("up", "down"):
            vals = scan_values(central, step, lo, hi, direction)
            n = len(vals)
            g_arr = vals if name == "grid factor" else np.full(n, GRID_INDIA)
            s_arr = vals if name != "grid factor" else np.full(n, SEC_H2DRI_0SCRAP)
            bi = best_per_draw(abatement_all_plans(g_arr, s_arr, np.tile(np.array(GEI_2026), (n, 1))), B, S)
            changed = np.nonzero(bi != ALL_PLANS.index(p0))[0]
            limit_txt = f"{(hi if direction == 'up' else lo):.6f}"
            if len(changed) == 0:
                flip_rows.append([case, f"{fB:.2f}", f"{fS:.2f}", name, unit, f"{central:.6f}", direction, f"{step:.6f}",
                                  limit_txt, plan_str(p0), f"{a0:.6f}", ru_txt, a_ru_txt, gap_txt, gap_pct_txt,
                                  "no change inside the search limit", "", "", "", ""])
                print(f"    {name:30s} {direction:4s}: no change of plan up to {limit_txt} {unit}")
                continue
            j = int(changed[0])
            v = float(vals[j])
            def exact_at(x):
                A1, C1, U1 = (tables(x, GEI_2026) if name == "grid factor" else tables(GRID_INDIA, GEI_2026, x))
                return solve(A1, C1, U1, B, S)
            p_new, a_new = exact_at(v)[0], exact_at(v)[1]
            assert p_new == ALL_PLANS[bi[j]] and p_new != p0, "exact solver does not confirm the flip point"
            v_before = central if j == 0 else float(vals[j - 1])
            assert exact_at(v_before)[0] == p0, "exact solver: the plan had already changed one grid step earlier"
            flip_rows.append([case, f"{fB:.2f}", f"{fS:.2f}", name, unit, f"{central:.6f}", direction, f"{step:.6f}",
                              limit_txt, plan_str(p0), f"{a0:.6f}", ru_txt, a_ru_txt, gap_txt, gap_pct_txt,
                              f"{v:.6f}", f"{v - central:.6f}", f"{100 * (v - central) / central:.6f}",
                              plan_str(p_new), f"{a_new:.6f}"])
            print(f"    {name:30s} {direction:4s}: plan changes at {v:.4f} {unit} ({v - central:+.4f}, "
                  f"{100 * (v - central) / central:+.2f} %) -> {plan_str(p_new)} ({a_new:.3f} Mt/yr there)")
with open(f"{OUT}/twolimit_flip_thresholds.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["limit_case", "f_B", "f_S", "input_moved", "unit", "central_value", "direction", "grid_step",
                 "search_limit", "best_plan_at_central", "abatement_at_central_Mt_per_yr", "runner_up_plan_at_central",
                 "runner_up_abatement_Mt_per_yr", "gap_best_minus_runner_up_Mt_per_yr", "gap_pct_of_best",
                 "first_grid_value_where_plan_changes", "change_from_central", "change_from_central_pct",
                 "new_best_plan", "new_best_plan_abatement_there_Mt_per_yr"])
    wr.writerows(flip_rows)

with open(f"{OUT}/software_versions_phase7_two_limit_selection_colab.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["package", "version"])
    wr.writerows([["python", platform.python_version()], ["numpy", np.__version__]])

print("\nFiles written in", OUT + "/: twolimit_plant_options.csv, twolimit_decisions.csv, twolimit_summary.csv, "
      "twolimit_uncertainty_plans.csv, twolimit_uncertainty_summary.csv, twolimit_flip_thresholds.csv")
print(f"Run time {time.time() - T0:.1f} s")
