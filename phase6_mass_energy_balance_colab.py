# -*- coding: utf-8 -*-
"""
phase6_mass_energy_balance_colab.py
===================================
Mass and energy balance of the hydrogen direct reduction + electric arc furnace route
(H2-DRI + EAF), per tonne of liquid steel (tLS), built ONLY from figures printed in
Vogl, Ahman & Nilsson (2018), J. Clean. Prod. 203, 736-745 (main text; the supplementary
file of that paper was not available). Where the paper gives no figure, the stream is
left out and the balance is reported as "partial". Nothing is taken from other sources
except standard atomic masses.

HOW TO RUN IN GOOGLE COLAB
  1. Run phase2_route_estimates_colab.py and phase3_route_model_colab.py first (same runtime).
     This script reads outputs_phase3/scenario_table.csv.
  2. Run:   !python phase6_mass_energy_balance_colab.py
Outputs (written to ./outputs_phase3/):
  mass_energy_balance.csv        one row per scrap share and stream (value, unit, source tag)
  balance_checks.csv             one row per check (expected, obtained, difference, status)
  mass_energy_balance_sweep.csv  fine grid, scrap share 0 to 1 in steps of 0.01
  software_versions_phase6_mass_energy_balance_colab.csv   python version of the run (no other package is imported)

Tags: [VERIFIED] read in the Vogl et al. (2018) text, "PDF PAGE n" = page of the PDF file
(PDF page 1 = journal page 736); [ASSUMPTION] study choice; [IMPLIED] not printed in the
paper, back-calculated here from printed figures; [COMPUTED] from the items above;
[CONSTANT] standard atomic mass.

Status rule for checks (fixed before the run, same for every check):
  pass         |relative difference| <= 1 %
  fail         |relative difference| >  1 %
  partial      the paper lacks at least one stream, so the balance cannot be closed
  not possible the paper gives no figure to compare with

DO NOT QUOTE THE NUMBER OF PASSED CHECKS AS A RESULT. Many checks are identities: they compare two
forms of the same formula, or a scaling with the row it was fitted to, and pass by construction (their
notes say so). The count of "pass" rows therefore says nothing about how well the balance closes. Quote
single checks, with their difference, and the "partial" and "not possible" rows.
"""

import os, csv

SEED = 12345            # no random numbers are used; kept for consistency with the other scripts
TOL = 0.01              # 1 % relative tolerance for "pass"                          [ASSUMPTION]
OUT = "outputs_phase3"
os.makedirs(OUT, exist_ok=True)

# =============================================================================
# PART 1 -- CONSTANTS
# =============================================================================
# Standard atomic masses (g/mol)                                                   [CONSTANT]
M_FE, M_O, M_H, M_C, M_CA = 55.845, 15.999, 1.008, 12.011, 40.078
M_FE2O3 = 2 * M_FE + 3 * M_O
M_H2, M_H2O, M_O2 = 2 * M_H, 2 * M_H + M_O, 2 * M_O
M_CO2, M_CAO = M_C + 2 * M_O, M_CA + M_O

# --- Vogl et al. (2018), compositions -------------------------------------------------
X_HEMATITE_PELLET = 0.95  # PDF PAGE 2: "The iron ore pellets considered contain 95% hematite (Fe2O3)
                          #   and 5% inert substances."                                        [VERIFIED]
X_FE_SCRAP = 0.95         # PDF PAGE 2: "Scrap charged to the EAF contains 95% iron and 5% inert
                          #   substances."                                                     [VERIFIED]
X_FE_STEEL = 1.00         # PDF PAGE 3: "The liquid steel product contains only iron, whereas all
                          #   other elements leave the EAF through the slag." Also: "all iron entering
                          #   the EAF leaves the process through the liquid steel product."   [VERIFIED]
# --- material flows, kg per tLS --------------------------------------------------------
PELLETS_0, H2_0 = 1504.0, 51.0
                          # PDF PAGE 4: "For continuous operation on 100% HBI, the production of one
                          #   tonne of steel requires 1504 kg of iron ore pellets. Without considering
                          #   hydrogen losses, 51 kg of hydrogen is needed per tonne of steel output."  [VERIFIED]
PELLETS_50, SCRAP_50, H2_50 = 738.0, 536.0, 25.0
                          # PDF PAGE 4: "If the EAF is charged with equal shares of HBI and scrap, 738 kg
                          #   of pellets, 536 kg of scrap and 25 kg of hydrogen are needed per tonne of
                          #   steel."                                                          [VERIFIED]
LIME = 50.0               # PDF PAGE 3: "a lime consumption in the EAF of 50 kg/tLS was assumed"; "the
                          #   flow to the EAF is considered to be pure CaO"                    [VERIFIED]
ELECTRODES = 2.0          # PDF PAGE 4: "graphite electrodes are consumed at a rate of 2 kg/tLS"  [VERIFIED]
ALLOYS = 11.0             # PDF PAGE 4: "Alloy use in the EAF was set at 11 kg per tonne of steel"
                          #   (cost model only; conflicts with "steel contains only iron")     [VERIFIED]
O2_VOGL = 411.0           # PDF PAGE 6: "For every tonne of liquid steel, 411 kg of oxygen are produced."
                          #   (text extraction prints "41 1")                                  [VERIFIED]
CARBON = None             # Mass of carbon added to the EAF: NOT in the main text. PDF PAGE 3 gives only
                          #   "half of the supplied carbon enters the steel, whereas the other half is
                          #   converted to CO2".                                               [NOT IN SOURCE]
METALLISATION = None      # Base-case value NOT in the main text. PDF PAGE 3: "only high metallisation
                          #   (>90%) is investigated"; "If the metallisation of HBI is lower than 94% this
                          #   additionally increases the SEC of the EAF".                      [NOT IN SOURCE]
LAMBDA_H2 = 1.5           # PDF PAGE 5: "For the base case a value of lambda = 1.5 was set." Feed ratio only;
                          #   surplus hydrogen is recycled, so it does not change the net need. [VERIFIED]
# --- energy ---------------------------------------------------------------------------
SEC_0 = 3.48              # PDF PAGE 4: "A specific energy consumption (SEC) of 3.48 MWh/t of liquid steel
                          #   was calculated for the base case."                               [VERIFIED]
EAF_SCRAP, EAF_DRI = 0.667, 0.753
                          # PDF PAGE 5: "producing steel in the EAF from scrap requires less energy
                          #   (0.667 MWh/tLS) than from pure DRI (0.753 MWh/tLS)"              [VERIFIED]
                          # PDF PAGE 3: "Both dependencies (HBI charge, metallisation) were implemented
                          #   as linear functions."                                            [VERIFIED]
ELY_SHARE = 2.0 / 3.0     # PDF PAGE 4: "The electrolyser consumes two thirds of the energy"   [VERIFIED]
ELY_EFF_LHV = 0.72        # PDF PAGE 3: "For the electrolyser an efficiency based on the lower heating
                          #   value (LHV) of 72% was assumed."                                 [VERIFIED]
DH_H2O_KJ_MOL = 242.0     # PDF PAGE 3, Eq. (1): "H2O(g) -> H2 + 1/2 O2, DHR = +242 kJ/mol". Used here as the
                          #   LHV of hydrogen per mole; the paper prints no LHV in MJ/kg.      [VERIFIED]
DH_RED_KJ_MOL = 99.5      # PDF PAGE 3, Eq. (2): "Fe2O3 + 3H2 -> 2Fe + 3H2O(g), DHR = +99.5 kJ/mol"  [VERIFIED]
PLANT_MT, ELY_MW, EAF_MW, PLANT_TWH, FREED_MW = 1.5, 411.0, 129.0, 5.22, 468.0
                          # PDF PAGE 8: "a plant capacity of 1.5 Mt/y that for optimal continuous operation
                          #   needs 129 MW of installed EAF capacity and 411 MW in electrolyser, resulting in
                          #   an electricity consumption of 5.22 TWh/y"                        [VERIFIED]
                          # PDF PAGE 9: "A switch from pure DRI to pure scrap operation frees up 468 MW, if
                          #   both the electrolyser and the shaft are shut down."              [VERIFIED]
HOURS = 8760.0            # continuous operation, all hours of the year. PDF PAGE 4 uses "8760 operating
                          #   hours" for the electrolyser lifetime; for the plant it is our reading. [ASSUMPTION]
GEI_BE_0, GEI_BE_25 = 532.0, 661.0
                          # PDF PAGE 5: "At pure HBI operation, the break-even GEI is 532 kgCO2/MWh." "At a
                          #   25% scrap charge in the EAF, the break-even GEI reaches 661 kgCO2/MWh." [VERIFIED]
NATIONAL = [("Germany", 29.5, 103.0), ("Sweden", 3.10, 10.8)]
                          # PDF PAGE 4: "German (29.5 Mt/y) or Swedish (3.10 Mt/y) primary steel production
                          #   ... would require about 103 TWh and 10.8 TWh of electricity"     [VERIFIED]
# --- CO2 ------------------------------------------------------------------------------
DIRECT_CO2_KG = 53.0      # PDF PAGE 6: "Emissions from carbon and lime use and consumption of the graphite
                          #   electrodes would result in emissions of 53 kg CO2 per tonne of steel"  [VERIFIED]
                          # PDF PAGE 3: "1 mol of CO2 is created for each mole of calcium oxide"; "Graphite
                          #   electrodes emit 1 mol of CO2 per mole of carbon that is consumed"  [VERIFIED]
# --- route model of Step III (phase3_route_model_colab.py; constants repeated, not changed) ----
REST_FUEL_GJ = 2.2        # PDF PAGE 4-5: "For pelletizing and rolling an additional 2.2 GJ/t of fuel as
REST_ELEC_MWH = 0.140     #   well as 140 kWh/t electricity are required irrespective of the production
                          #   route (Worrell et al., 2008)."                                   [VERIFIED]
EF_NATGAS = 0.0561        # t CO2/GJ. IPCC 2006 Guidelines, Vol. 2, Ch. 1, Table 1.4
                          #   (PDF PAGE 24): "Natural Gas 15.3 1 56 100 54 300 58 300" -> 56,100 kg CO2/TJ.
                          #   Same value as in phase3_route_model_colab.py.                    [VERIFIED]
SCRAP_SHARES = [0.0, 0.50, 0.95, 1.00]   # the three shares of the study; 1.00 added to show the end point of the
                                         # straight line                                       [ASSUMPTION]

# =============================================================================
# PART 2 -- DERIVED QUANTITIES
# =============================================================================
FE_IN_FE2O3 = 2 * M_FE / M_FE2O3                    # 0.6994 (pure hematite)
FE_IN_PELLET = X_HEMATITE_PELLET * FE_IN_FE2O3      # iron content of Vogl's pellets            [COMPUTED]
O_IN_PELLET = X_HEMATITE_PELLET * (1 - FE_IN_FE2O3)
O_PER_FE_FEO = M_O / M_FE                           # kg O per kg Fe bound as FeO
LHV_KWH_KG = DH_H2O_KJ_MOL / M_H2 / 3.6             # kWh per kg H2 (LHV from Vogl Eq. 1)       [COMPUTED]
ELY_KWH_KG = LHV_KWH_KG / ELY_EFF_LHV               # electrolyser electricity per kg H2        [COMPUTED]
H2_STOICH_PER_T_FE = 1000 * 3 * M_H2 / (2 * M_FE)   # kg H2 per t Fe, Fe2O3 + 3H2 -> 2Fe + 3H2O [COMPUTED]

def hbi_per_kg_pellet(m):
    """kg HBI per kg pellet at metallisation m (iron + inert + oxygen left in FeO)."""
    return FE_IN_PELLET * (1 + (1 - m) * O_PER_FE_FEO) + (1 - X_HEMATITE_PELLET)

def h2_need(fe_ore_kg, m):
    """kg H2 to bring fe_ore_kg of iron from Fe2O3 to metallisation m (rest stays as FeO)."""
    return fe_ore_kg / M_FE * (1.5 * m + 0.5 * (1 - m)) * M_H2

# Metallisation is not printed. Three back-calculations from printed figures:              [IMPLIED]
#  (a) 50% row: "equal shares of HBI and scrap" read as equal MASS -> HBI = 536 kg from 738 kg pellets
M_IMPL_MASS = 1 - (SCRAP_50 / PELLETS_50 - hbi_per_kg_pellet(1.0)) / (FE_IN_PELLET * O_PER_FE_FEO)
#  (b) oxygen: 411 kg O2 leaves; the rest of the ore oxygen stays in FeO
M_IMPL_O2 = 1 - (PELLETS_0 * O_IN_PELLET - O2_VOGL) / O_PER_FE_FEO / (PELLETS_0 * FE_IN_PELLET)
#  (c) hydrogen: 51 kg H2 for the iron in 1504 kg pellets
M_IMPL_H2 = (H2_0 / M_H2 / (PELLETS_0 * FE_IN_PELLET / M_FE) - 0.5) / 1.0
M_USED = M_IMPL_MASS      # used for the iron-balance scaling below; (a) is the only one that fixes HBI mass
H_HBI = hbi_per_kg_pellet(M_USED)

# --- two ways to scale the streams with the scrap share s ------------------------------
# P: proportional to the ore share (1 - s). This is what the straight electricity line of the route model
#    implies. It is kept here for comparison only: the plant-choice model (phase 7, phase 8) does NOT use
#    1072 x s for scrap; it reads the scrap mass of the iron-balance scaling M below.     [ASSUMPTION]
def pellets_P(s): return PELLETS_0 * (1 - s)
def scrap_P(s):   return SCRAP_50 / 0.5 * s
def h2_P(s):      return H2_0 * (1 - s)
# M: s = scrap mass / (scrap mass + HBI mass), iron balance closed on 1000 kg Fe in steel, compositions
#    as printed, metallisation M_USED. Reading "share" as a mass share is our reading.     [ASSUMPTION]
def charge_M(s):
    total = 1000.0 / ((1 - s) * FE_IN_PELLET / H_HBI + s * X_FE_SCRAP)   # kg HBI + scrap per tLS
    return (1 - s) * total / H_HBI, s * total                             # pellets, scrap
def pellets_M(s): return charge_M(s)[0]
def scrap_M(s):   return charge_M(s)[1]
def h2_M(s):      return H2_0 * pellets_M(s) / pellets_M(0.0)            # hydrogen follows the ore
def ore_ratio_M(s): return pellets_M(s) / pellets_M(0.0)

# --- energy ---------------------------------------------------------------------------
def sec_line(s):          # straight line of the route model (phase3_route_model_colab.py)
    return SEC_0 - s * (SEC_0 - EAF_SCRAP)
def eaf(s):               # linear in the HBI share of the charge (Vogl, PDF PAGE 3)        [VERIFIED form]
    return EAF_DRI - s * (EAF_DRI - EAF_SCRAP)
ELY_0_FROM_H2 = H2_0 * ELY_KWH_KG / 1000              # MWh/tLS, from 51 kg and 72 % LHV    [COMPUTED]
ELY_0_FROM_MW = ELY_MW * HOURS / (PLANT_MT * 1e6)     # MWh/tLS, from 411 MW                [COMPUTED]
ELY_0_FROM_SHARE = ELY_SHARE * SEC_0                  # MWh/tLS, from "two thirds"          [COMPUTED]
UPSTREAM_0 = SEC_0 - EAF_DRI                          # everything before the EAF at 0 % scrap
OTHER_0 = UPSTREAM_0 - ELY_0_FROM_H2                  # ore heating, shaft, hydrogen heating: the paper prints
                                                      # no number (Fig. 2 only), so this is a remainder  [IMPLIED]
def sec_bottom_up(s, ratio):
    """EAF (linear, Vogl) + everything before the EAF scaled with the ore ratio.          [ASSUMPTION]
    Vogl, PDF PAGE 5: "less iron ore is required per output, which leads to lower energy consumption
    in all process steps before the EAF"."""
    return eaf(s) + UPSTREAM_0 * ratio

def intensity(sec_value, g):    # t CO2 per tLS, same form as phase3_route_model_colab.py
    return DIRECT_CO2_KG / 1000 + sec_value * g + REST_FUEL_GJ * EF_NATGAS + REST_ELEC_MWH * g

# =============================================================================
# PART 3 -- STREAM TABLE
# =============================================================================
rows = []
def add(s, group, stream, value, unit, tag):
    rows.append([f"{s:.2f}", group, stream, "" if value is None else f"{value:.6f}", unit, tag])

for s in SCRAP_SHARES:
    pM, scM, hM = pellets_M(s), scrap_M(s), h2_M(s)
    pP, scP, hP = pellets_P(s), scrap_P(s), h2_P(s)
    fe_ore = pM * FE_IN_PELLET
    o_removed = hM * M_O / M_H2                      # each H2 removes one O atom
    o_left_feo = pM * O_IN_PELLET - o_removed        # oxygen still bound in the HBI
    inert = pM * (1 - X_HEMATITE_PELLET) + scM * (1 - X_FE_SCRAP)
    # inputs
    add(s, "input", "iron ore pellets (iron-balance scaling M)", pM, "kg/tLS", "COMPUTED from VERIFIED compositions; s=0, 0.5 printed as 1504, 738")
    add(s, "input", "scrap (iron-balance scaling M)", scM, "kg/tLS", "COMPUTED; s=0.5 printed as 536")
    add(s, "input", "iron ore pellets (proportional scaling P)", pP, "kg/tLS", "ASSUMPTION 1504 x (1 - s)")
    add(s, "input", "scrap (proportional scaling P)", scP, "kg/tLS", "ASSUMPTION 1072 x s")
    add(s, "input", "lime (CaO)", LIME, "kg/tLS", "VERIFIED p3; per tLS at every scrap share")
    add(s, "input", "graphite electrodes", ELECTRODES, "kg/tLS", "VERIFIED p4; per tLS at every scrap share")
    add(s, "input", "alloys", ALLOYS, "kg/tLS", "VERIFIED p4 (cost model); not in Vogl's mass balance")
    add(s, "input", "carbon added to EAF", None, "kg/tLS", "NOT IN SOURCE (main text)")
    # internal loop
    add(s, "internal loop", "hydrogen consumed (scaling M)", hM, "kg/tLS", "COMPUTED 51 x ore ratio; s=0, 0.5 printed as 51, 25")
    add(s, "internal loop", "hydrogen consumed (scaling P)", hP, "kg/tLS", "ASSUMPTION 51 x (1 - s)")
    add(s, "internal loop", "water from reduction, returned to electrolyser", hM * M_H2O / M_H2, "kg/tLS", "COMPUTED stoichiometry; no water crosses the boundary (p2)")
    # outputs
    add(s, "output", "liquid steel", 1000.0, "kg/tLS", "basis; contains only iron (VERIFIED p3)")
    add(s, "output", "oxygen from electrolyser", o_removed, "kg/tLS", "COMPUTED from hydrogen; s=0 printed as 411")
    add(s, "output", "slag: inert from pellets and scrap + lime (lower bound)", inert + LIME, "kg/tLS", "COMPUTED; slag total NOT IN SOURCE -> not closed")
    add(s, "output", "oxygen left in FeO of the HBI (fate in EAF not stated)", o_left_feo, "kg/tLS", "COMPUTED remainder; depends on hydrogen figure")
    add(s, "output", "CO2 from lime (calcination, 1 mol per mol CaO)", LIME * M_CO2 / M_CAO, "kg/tLS", "COMPUTED from VERIFIED p3")
    add(s, "output", "CO2 from electrodes", ELECTRODES * M_CO2 / M_C, "kg/tLS", "COMPUTED from VERIFIED p3-4")
    add(s, "output", "CO2 from added carbon", None, "kg/tLS", "NOT IN SOURCE (carbon mass missing)")
    # iron
    add(s, "iron", "iron in pellets (M)", fe_ore, "kg Fe/tLS", "COMPUTED 95 % hematite")
    add(s, "iron", "iron in scrap (M)", scM * X_FE_SCRAP, "kg Fe/tLS", "COMPUTED 95 % iron")
    # electricity
    ely = hM * ELY_KWH_KG / 1000
    bu = sec_bottom_up(s, ore_ratio_M(s))
    add(s, "electricity", "electrolyser (hydrogen x LHV / 0.72)", ely, "MWh/tLS", "COMPUTED from VERIFIED p3-4")
    add(s, "electricity", "EAF (linear 0.753 -> 0.667)", eaf(s), "MWh/tLS", "VERIFIED p3, p5")
    add(s, "electricity", "ore heating, shaft, hydrogen heating (remainder x ore ratio)", bu - eaf(s) - ely, "MWh/tLS", "IMPLIED remainder; no printed number")
    add(s, "electricity", "total, bottom-up (M)", bu, "MWh/tLS", "COMPUTED")
    add(s, "electricity", "total, straight line of the route model", sec_line(s), "MWh/tLS", "phase3_route_model_colab.py")
    add(s, "electricity", "bottom-up minus straight line", bu - sec_line(s), "MWh/tLS", "COMPUTED")
    add(s, "heat (information)", "reaction heat of reduction, Eq. 2, full reduction of the ore iron", fe_ore * 1000 / M_FE / 2 * DH_RED_KJ_MOL / 3.6e6, "MWh/tLS", "COMPUTED from VERIFIED p3; part of the remainder, not added")

with open(f"{OUT}/mass_energy_balance.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["scrap_share", "group", "stream", "value", "unit", "source_tag"])
    wr.writerows(rows)

# =============================================================================
# PART 4 -- CHECKS
# =============================================================================
checks = []
def chk(name, expected, obtained, unit, note, status=None):
    diff = obtained - expected
    rel = diff / expected if expected else float("nan")
    if status is None:
        status = "pass" if abs(rel) <= TOL else "fail"
    checks.append([name, f"{expected:.6f}", f"{obtained:.6f}", f"{diff:.6f}", f"{100 * rel:.4f}", unit, status, note])
def na(name, status, note):
    checks.append([name, "", "", "", "", "", status, note])

# --- iron balance on the printed rows ---
chk("iron balance, 0% scrap (printed 1504 kg pellets)", 1000.0, PELLETS_0 * FE_IN_PELLET, "kg Fe/tLS",
    "expected = iron in 1 t steel (steel = iron only); obtained = iron in pellets at 95 % hematite")
chk("iron balance, 50% scrap (printed 738 kg pellets + 536 kg scrap)", 1000.0,
    PELLETS_50 * FE_IN_PELLET + SCRAP_50 * X_FE_SCRAP, "kg Fe/tLS", "pellets 95 % hematite, scrap 95 % iron")
chk("iron content of pellets: implied by 1504 kg/tLS vs stated composition", FE_IN_PELLET, 1000.0 / PELLETS_0,
    "kg Fe/kg pellet", "stated: 95 % of pure hematite (69.94 % Fe) = 66.45 % Fe; implied below 69.94 %, so physically possible")
# --- scaling residuals at the printed 50% row ---
chk("scaling P vs printed 50% row: pellets", PELLETS_50, pellets_P(0.5), "kg/tLS", "P = 1504 x (1 - s)")
chk("scaling P vs printed 50% row: hydrogen", H2_50, h2_P(0.5), "kg/tLS", "P = 51 x (1 - s); printed value has 2 digits (24.5-25.5)")
chk("scaling P vs printed 50% row: scrap", SCRAP_50, scrap_P(0.5), "kg/tLS", "zero by construction (1072 x s is fitted to this row)")
chk("scaling M vs printed 0% row: pellets", PELLETS_0, pellets_M(0.0), "kg/tLS", "M = iron balance; metallisation implied from the 50% row")
chk("scaling M vs printed 50% row: pellets", PELLETS_50, pellets_M(0.5), "kg/tLS", "M")
chk("scaling M vs printed 50% row: scrap", SCRAP_50, scrap_M(0.5), "kg/tLS", "M")
chk("scaling M vs printed 50% row: hydrogen", H2_50, h2_M(0.5), "kg/tLS", "M; printed value has 2 digits")
for s in SCRAP_SHARES:
    chk(f"iron balance of scaling P at {int(round(100 * s))}% scrap", 1000.0,
        pellets_P(s) * FE_IN_PELLET + scrap_P(s) * X_FE_SCRAP, "kg Fe/tLS",
        "P = 1504 x (1 - s) pellets and 1072 x s scrap (comparison only; the plant-choice model now uses the iron-balance scaling M)")
chk("scrap per tLS at 95% scrap: 1.072 x s vs iron-balance scaling M", scrap_M(0.95), scrap_P(0.95), "kg/tLS",
    "expected = M; obtained = 1072 x s")
chk("scrap per tLS at 100% scrap: 1.072 x s vs 1000 / 0.95", 1000.0 / X_FE_SCRAP, scrap_P(1.0), "kg/tLS",
    "pure scrap at 95 % iron needs 1052.6 kg; 1072 kg carries more iron than the steel holds")
# --- hydrogen ---
fe0 = PELLETS_0 * FE_IN_PELLET
chk("hydrogen: stoichiometric need per t iron (Fe2O3 + 3H2) vs printed 51 kg/tLS", H2_STOICH_PER_T_FE, H2_0, "kg H2",
    "ratio obtained/expected = share of full reduction; below 1 because the HBI is not fully metallised")
chk("hydrogen: need for the iron in 1504 kg pellets at implied metallisation (a) vs printed 51", h2_need(fe0, M_IMPL_MASS), H2_0, "kg H2/tLS",
    f"metallisation (a) = {M_IMPL_MASS:.4f} from the 50% row masses [IMPLIED]")
chk("hydrogen: implied by 411 kg oxygen vs printed 51", O2_VOGL * M_H2 / M_O, H2_0, "kg H2/tLS",
    "expected = 411 x 2.016 / 15.999; one H2 per O atom")
chk("oxygen from electrolyser: from 51 kg hydrogen vs printed 411", O2_VOGL, H2_0 * M_O / M_H2, "kg O2/tLS", "same comparison, other direction")
chk("hydrogen ratio 50% / 0% vs pellet ratio 738 / 1504", PELLETS_50 / PELLETS_0, H2_50 / H2_0, "-", "hydrogen follows the ore")
chk("metallisation implied by oxygen (b) vs implied by masses (a)", M_IMPL_MASS, M_IMPL_O2, "-", "both IMPLIED; the paper prints no base-case value", )
chk("metallisation implied by hydrogen (c) vs implied by masses (a)", M_IMPL_MASS, M_IMPL_H2, "-", "both IMPLIED; (c) rests on the 2-digit figure 51")
# --- energy ---
chk("total: 5.22 TWh/y / 1.5 Mt/y vs printed 3.48", SEC_0, PLANT_TWH / PLANT_MT, "MWh/tLS", "plant example, PDF page 8")
for name, mt, twh in NATIONAL:
    chk(f"total: {name} {twh} TWh / {mt} Mt vs printed 3.48", SEC_0, twh / mt, "MWh/tLS", "PDF page 4; 'about'")
chk("EAF: 129 MW x 8760 h / 1.5 Mt vs printed 0.753", EAF_DRI, EAF_MW * HOURS / (PLANT_MT * 1e6), "MWh/tLS", "8760 h/y is an assumption")
chk("electrolyser: 51 kg x LHV / 0.72 vs 'two thirds' of 3.48", ELY_0_FROM_SHARE, ELY_0_FROM_H2, "MWh/tLS",
    f"LHV {LHV_KWH_KG:.3f} kWh/kg from Eq. 1; {ELY_KWH_KG:.3f} kWh/kg electricity; 'two thirds' is a rounded phrase")
chk("electrolyser: 51 kg x LHV / 0.72 vs 411 MW x 8760 h / 1.5 Mt", ELY_0_FROM_MW, ELY_0_FROM_H2, "MWh/tLS", "8760 h/y is an assumption")
chk("electrolyser: 411 MW x 8760 h / 1.5 Mt vs 'two thirds' of 3.48", ELY_0_FROM_SHARE, ELY_0_FROM_MW, "MWh/tLS", "")
chk("sum of parts at 0% scrap vs 3.48 (electrolyser from 51 kg + EAF 0.753; ore heating and shaft missing)",
    SEC_0, ELY_0_FROM_H2 + EAF_DRI, "MWh/tLS",
    "the difference is the energy for ore heating, shaft and hydrogen heating; the paper gives it only in Fig. 2 -> cannot close",
    status="partial")
chk("freed power, pure DRI -> pure scrap: 468 MW vs straight-line slope 3.48 - 0.667", SEC_0 - EAF_SCRAP,
    FREED_MW * HOURS / (PLANT_MT * 1e6), "MWh/tLS", "PDF page 9; 8760 h/y assumed")
chk("freed power: 468 MW vs 3.48 - 0.753 (EAF power kept at the DRI value)", SEC_0 - EAF_DRI,
    FREED_MW * HOURS / (PLANT_MT * 1e6), "MWh/tLS", "the 468 MW equals electrolyser + other upstream only")
SEC_25_IMPLIED = SEC_0 * GEI_BE_0 / GEI_BE_25
chk("straight line at 25% scrap vs value implied by Vogl's break-even pair 532 / 661", SEC_25_IMPLIED, sec_line(0.25), "MWh/tLS",
    "expected = 3.48 x 532 / 661; valid only if Vogl's non-electric emissions are the same at 0 % and 25 % scrap [ASSUMPTION]")
chk("bottom-up (M) at 25% scrap vs value implied by Vogl's break-even pair 532 / 661", SEC_25_IMPLIED,
    sec_bottom_up(0.25, ore_ratio_M(0.25)), "MWh/tLS", "same assumption")
for s in SCRAP_SHARES:
    chk(f"electricity at {int(round(100 * s))}% scrap: bottom-up with scaling P vs straight line", sec_line(s),
        sec_bottom_up(s, 1 - s), "MWh/tLS", "identical by algebra: both are linear in s with the same end points")
for s in SCRAP_SHARES:
    chk(f"electricity at {int(round(100 * s))}% scrap: bottom-up with scaling M vs straight line", sec_line(s),
        sec_bottom_up(s, ore_ratio_M(s)), "MWh/tLS", "expected = straight line of the existing model; obtained = bottom-up")
chk("electricity at 50% scrap: bottom-up with printed ore ratio 738/1504 vs straight line", sec_line(0.5),
    sec_bottom_up(0.5, PELLETS_50 / PELLETS_0), "MWh/tLS", "uses only printed figures and the ore-ratio assumption")
na("energy for ore heating, shaft and hydrogen heating, computed from heat capacities", "not possible",
   "heat capacities, condenser recovery and the Fig. 2 values are not in the main text")
# --- CO2 ---
co2_lime, co2_el = LIME * M_CO2 / M_CAO, ELECTRODES * M_CO2 / M_C
chk("direct CO2: lime + electrodes vs printed 53 kg/tLS (added carbon missing)", DIRECT_CO2_KG, co2_lime + co2_el, "kg CO2/tLS",
    "carbon mass is not in the main text; the difference equals the CO2 from added carbon if the paper is consistent", status="partial")
C_IMPLIED = (DIRECT_CO2_KG - co2_lime - co2_el) * M_C / M_CO2 * 2
na("carbon added to the EAF", "not possible",
   f"not printed; implied {C_IMPLIED:.3f} kg/tLS if half becomes CO2 (PDF page 3) [IMPLIED, not a source value]")
# --- overall mass ---
m_in = PELLETS_0 + LIME + ELECTRODES
m_out = 1000.0 + O2_VOGL + PELLETS_0 * (1 - X_HEMATITE_PELLET) + LIME + ELECTRODES
chk("overall mass at 0% scrap, printed figures: inputs (pellets, lime, electrodes) vs outputs (steel, 411 kg O2, inert + lime slag, electrode carbon)",
    m_in, m_out, "kg/tLS", "not counted: added carbon, alloys, oxygen left in FeO, air oxygen in CO2; slag total not in source", status="partial")
m_in50 = PELLETS_50 + SCRAP_50 + LIME + ELECTRODES
m_out50 = 1000.0 + H2_50 * M_O / M_H2 + PELLETS_50 * 0.05 + SCRAP_50 * 0.05 + LIME + ELECTRODES
chk("overall mass at 50% scrap, printed figures: inputs vs outputs (oxygen computed from 25 kg hydrogen)",
    m_in50, m_out50, "kg/tLS", "same omissions; no printed oxygen figure at 50 % scrap", status="partial")

# --- consistency with phase3 scenario_table.csv ---
SCEN = f"{OUT}/scenario_table.csv"
if not os.path.exists(SCEN):
    raise FileNotFoundError("Run phase3_route_model_colab.py first: " + SCEN + " not found.")
with open(SCEN, encoding="utf-8") as f:
    scen = list(csv.DictReader(f))
GRID_INDIA = None
for r in scen:
    s, g = float(r["scrap_share"]), float(r["grid_gCO2_per_kWh"]) / 1000
    if r["electricity"] == "India grid":
        GRID_INDIA = g
    tag = f"{int(round(100 * s))}% scrap, {r['electricity']}"
    chk(f"phase 3 SEC reproduced by straight line, {tag}", float(r["SEC_MWh_per_t"]), sec_line(s), "MWh/tLS", "scenario_table.csv")
    chk(f"phase 3 route intensity reproduced with straight-line electricity, {tag}", float(r["I_new_tCO2_per_t"]),
        intensity(sec_line(s), g), "t CO2/tLS", "scenario_table.csv")
    chk(f"route intensity with bottom-up (M) electricity vs phase 3 value, {tag}", float(r["I_new_tCO2_per_t"]),
        intensity(sec_bottom_up(s, ore_ratio_M(s)), g), "t CO2/tLS", "expected = existing model; obtained = bottom-up; both kept")
assert GRID_INDIA is not None, "scenario_table.csv has no 'India grid' row"

with open(f"{OUT}/balance_checks.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["check", "expected", "obtained", "difference", "relative_difference_pct", "unit", "status", "note"])
    wr.writerows(checks)

# =============================================================================
# PART 5 -- FINE GRID
# =============================================================================
with open(f"{OUT}/mass_energy_balance_sweep.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["scrap_share", "pellets_P_kg", "pellets_M_kg", "scrap_P_kg", "scrap_M_kg", "hydrogen_P_kg", "hydrogen_M_kg",
                 "iron_in_P_kg", "iron_in_M_kg", "EAF_MWh", "SEC_straight_line_MWh", "SEC_bottom_up_M_MWh",
                 "SEC_difference_MWh", "SEC_difference_pct", "I_grid_straight_line_tCO2", "I_grid_bottom_up_M_tCO2",
                 "I_grid_difference_tCO2"])
    for k in range(101):
        s = k / 100
        bu, ln = sec_bottom_up(s, ore_ratio_M(s)), sec_line(s)
        vals = [pellets_P(s), pellets_M(s), scrap_P(s), scrap_M(s), h2_P(s), h2_M(s),
                pellets_P(s) * FE_IN_PELLET + scrap_P(s) * X_FE_SCRAP, pellets_M(s) * FE_IN_PELLET + scrap_M(s) * X_FE_SCRAP,
                eaf(s), ln, bu, bu - ln, 100 * (bu - ln) / ln, intensity(ln, GRID_INDIA), intensity(bu, GRID_INDIA),
                intensity(bu, GRID_INDIA) - intensity(ln, GRID_INDIA)]
        wr.writerow([f"{s:.2f}"] + [f"{v:.6f}" for v in vals])

# =============================================================================
# PART 6 -- ASSERTS ON REAL PROPERTIES, AND PRINTS
# =============================================================================
for k in range(101):
    s = k / 100
    assert abs(scrap_M(s) - s * (scrap_M(s) + pellets_M(s) * H_HBI)) < 1e-9, "s must be the scrap mass share of the charge"
    assert pellets_M(s) >= -1e-12 and h2_M(s) >= -1e-12
assert abs(M_FE2O3 - 159.687) < 1e-9 and 0.69 < FE_IN_FE2O3 < 0.70
assert 0.9 < M_USED < 1.0, "implied metallisation outside the range the paper studies (>90 %)"

print("=" * 100)
print("Derived constants")
print(f"  iron in pure Fe2O3 {100 * FE_IN_FE2O3:.3f} %; in Vogl's pellets {100 * FE_IN_PELLET:.3f} %; implied by 1504 kg {100 * 1000 / PELLETS_0:.3f} %")
print(f"  H2 stoichiometric {H2_STOICH_PER_T_FE:.3f} kg/t Fe; printed 51 -> ratio {H2_0 / H2_STOICH_PER_T_FE:.4f}")
print(f"  LHV from Eq. 1 {LHV_KWH_KG:.3f} kWh/kg; electrolyser {ELY_KWH_KG:.3f} kWh/kg H2")
print(f"  electrolyser at 0% scrap: from 51 kg {ELY_0_FROM_H2:.4f}; from 411 MW {ELY_0_FROM_MW:.4f}; two thirds {ELY_0_FROM_SHARE:.4f} MWh/tLS")
print(f"  remainder (ore heating, shaft, H2 heating) {OTHER_0:.4f} MWh/tLS [IMPLIED]")
print(f"  implied metallisation: masses {M_IMPL_MASS:.4f}; oxygen {M_IMPL_O2:.4f}; hydrogen {M_IMPL_H2:.4f}")
print(f"  CO2 lime {co2_lime:.3f} + electrodes {co2_el:.3f} = {co2_lime + co2_el:.3f} kg vs 53; implied carbon {C_IMPLIED:.3f} kg/tLS")
print(f"  SEC at 25% scrap implied by 532/661: {SEC_25_IMPLIED:.4f}; straight line {sec_line(0.25):.4f}; bottom-up M {sec_bottom_up(0.25, ore_ratio_M(0.25)):.4f}")
print("\nStreams and electricity per tLS")
print(f"{'s':>5} {'pel_M':>8} {'pel_P':>8} {'scr_M':>8} {'scr_P':>8} {'H2_M':>6} {'EAF':>6} {'line':>7} {'b-up M':>7} {'diff':>8} {'%':>6} {'I line':>7} {'I b-up':>7}")
for s in SCRAP_SHARES:
    bu, ln = sec_bottom_up(s, ore_ratio_M(s)), sec_line(s)
    print(f"{s:5.2f} {pellets_M(s):8.1f} {pellets_P(s):8.1f} {scrap_M(s):8.1f} {scrap_P(s):8.1f} {h2_M(s):6.2f} {eaf(s):6.3f} "
          f"{ln:7.4f} {bu:7.4f} {bu - ln:+8.4f} {100 * (bu - ln) / ln:+6.2f} {intensity(ln, GRID_INDIA):7.4f} {intensity(bu, GRID_INDIA):7.4f}")
print("\nChecks")
for c in checks:
    print(f"  [{c[6]:>12}] {c[0][:92]:92s} exp {c[1]:>12} got {c[2]:>12} rel% {c[4]:>9}")
cnt = {}
for c in checks:
    cnt[c[6]] = cnt.get(c[6], 0) + 1
print("\nStatus counts:", cnt)
print("NOTE: do not quote the number of passed checks as a result; many checks are identities (see the docstring).")

import platform
with open(f"{OUT}/software_versions_phase6_mass_energy_balance_colab.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["package", "version"])
    wr.writerow(["python", platform.python_version()])
print("Phase 6 done. Files in", OUT)
