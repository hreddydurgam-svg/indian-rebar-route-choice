# -*- coding: utf-8 -*-
"""
phase3_route_model_colab.py
===========================
Step III - Route model with plant baselines: hydrogen direct reduction (H2-DRI) + electric arc
furnace (EAF) with a variable share of scrap, under two electricity cases (renewable, Indian grid),
compared with EACH PLANT'S OWN baseline intensity (no national average is used).

HOW TO RUN IN GOOGLE COLAB
  1. Run phase2_route_estimates_colab.py first (same runtime). It writes
     outputs_phase2/route_estimates.csv and outputs_phase2/epd_inputs.csv, which this script reads
     for ONE side-by-side table (route_vs_epd_evidence.csv). No EPD value enters any calculation here.
  2. Run:   !python phase3_route_model_colab.py
Outputs are written to ./outputs_phase3/
  plant_inputs.csv        plant, registration, output and GEI (2026 revised draft; 2025 draft)
  scenario_table.csv      scrap_share, electricity, grid_gCO2_per_kWh, SEC_MWh_per_t, I_new_tCO2_per_t
  plant_breakeven.csv     per plant: scrap share where I(s, India grid) = plant GEI; break-even grid factor at s = 0, 0.50, 0.95
  route_vs_epd_evidence.csv  side by side, NOT combined: route model intensity at s = 0, 0.50, 0.95 for both
                          electricity cases; and the two Step II route groups (input range and estimates)
  route_sensitivity.csv   main; grid 0.710 (unverified); without the pelletizing and rolling terms
  sweep_scrap.csv, sweep_grid.csv
  software_versions_phase3_route_model_colab.csv   python, numpy, matplotlib versions of the run
  fig_scrap_sweep.tiff, fig_grid_sweep.tiff   (1000 dpi + PDF, LZW)

Model (per tonne of liquid steel, t CO2 / t):
  SEC(s)       = SEC_0 - s*(SEC_0 - SEC_100)                      (linear in scrap share s)   [ASSUMPTION]
  I_new(s, g)  = DIRECT + SEC(s)*g + REST_FUEL*EF_NG + REST_ELEC*g
Per plant i (GEI_i = the plant's own gazette intensity, t CO2e / t of equivalent product):
  change_i     = I_new / GEI_i - 1
  break-even scrap share  : s where I_new(s, India grid) = GEI_i     (analytic inverse)
  break-even grid factor  : g where I_new(s, g) = GEI_i              (for s = 0, 0.50, 0.95)

NO EPD VALUE IS PAIRED WITH A PLANT GEI
  No plant's EPD value is scaled by I_new / GEI_i to get a "future rebar value".
  Reason: an EPD (cradle to gate, product level, CO2e) and a GEI (plant level, per tonne of
  equivalent product) come from different accounting systems, and their ratio is below 1 for two plants, so
  scaling one by the other gives values below the route intensity itself. This script therefore gives NO
  rebar value after conversion. In its place, route_vs_epd_evidence.csv shows the route model and the Step II
  EPD evidence side by side, without a ratio and without scaling, only to see whether the direction of the
  model with scrap share agrees with declared values.

LIMITS TO STATE IN THE PAPER
  * GEI (per tonne of equivalent product, plant boundary of the gazette) and I_new (per tonne of liquid
    steel, boundary of Vogl et al.) do not have the same unit basis or boundary. The comparison is indicative.
  * The route model (t CO2 per t steel; on-site and electricity only) and the EPD values (kg CO2e per kg rebar;
    cradle to gate) are different quantities. The side-by-side table is a direction check, not a validation.
  * Casting energy is not counted (Vogl et al. give no separate figure). Grid losses are not added.

Tags: [VERIFIED] read in the named source text (page marker and quote); [ASSUMPTION] study choice to state
in the paper; [UNVERIFIED] not checked against a primary source.
"""

# =============================================================================
# PART 0 -- COLAB DEPENDENCIES (runs only inside Colab)
# =============================================================================
import sys, subprocess, glob, os

IN_COLAB = "google.colab" in sys.modules or os.path.isdir("/content")
if IN_COLAB:
    for c in [
        "echo ttf-mscorefonts-installer msttcorefonts/accept-mscorefonts-eula select true | sudo debconf-set-selections",
        "sudo apt-get -qq update -y",
        "sudo DEBIAN_FRONTEND=noninteractive apt-get -qq install -y ttf-mscorefonts-installer",
        "rm -rf ~/.cache/matplotlib",
        f"{sys.executable} -m pip install -q matplotlib pillow numpy",
    ]:
        subprocess.run(c, shell=True, check=False)

# =============================================================================
# IMPORTS
# =============================================================================
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

for f in glob.glob("/usr/share/fonts/truetype/msttcorefonts/*.ttf"):
    try:
        fm.fontManager.addfont(f)
    except Exception:
        pass
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]

OUT = "outputs_phase3"
os.makedirs(OUT, exist_ok=True)
SEED = 12345          # not used: this script has no random step (kept for the common script layout)

def z0(v, nd=6):      # format without printing "-0.000000"
    v = round(float(v), nd)
    return f"{(0.0 if v == 0 else v):.{nd}f}"

# =============================================================================
# PART 1 -- PARAMETERS
# =============================================================================
# Vogl, Ahman & Nilsson (2018), J. Clean. Prod. 203:736-745 (1_s2_0_S0959652618326301_main_pdf.txt)
SEC_H2DRI_0SCRAP = 3.48   # MWh/t LS, H2-DRI + EAF, 100% HBI (0% scrap), base case. PDF PAGE 4: "A specific energy
                          # consumption (SEC) of 3.48 MWh/t of liquid steel was calculated for the base case."
                          # incl. electrolyser, shaft, EAF; excl. pelletizing, casting, rolling            [VERIFIED]
SEC_EAF_100SCRAP = 0.667  # MWh/t LS, EAF on 100% scrap. PDF PAGE 3: "A SEC of 0.667 MWh per tonne liquid steel
                          # (tLS) was used for the operation on pure scrap feed."                           [VERIFIED]
DIRECT_CO2       = 0.053  # t CO2/t LS. PDF PAGE 6: "Emissions from carbon and lime use and consumption of the
                          # graphite electrodes would result in emissions of 53 kg CO2 per tonne of steel"  [VERIFIED]
REST_FUEL_GJ     = 2.2    # GJ/t fuel and
REST_ELEC_MWH    = 0.140  # MWh/t electricity. PDF PAGE 4-5: "For pelletizing and rolling an additional 2.2 GJ/t of
                          # fuel as well as 140 kWh/t electricity are required irrespective of the production
                          # route (Worrell et al., 2008)."                                                  [VERIFIED]
                          # NOTE: Vogl gives no separate figure for casting -> casting energy not counted.

# IPCC 2006 Guidelines, Vol. 2, Ch. 1, Table 1.4 (PDF PAGE 24):
# "Natural Gas 15.3 1 56 100 54 300 58 300" -> 56,100 kg CO2/TJ (range 54,300-58,300)
EF_NATGAS        = 0.0561 # t CO2/GJ                                                              [VERIFIED]
                          # The 2.2 GJ fuel is natural gas: Vogl does not name the fuel           [ASSUMPTION]

# Central Electricity Authority, CO2 Baseline Database for the Indian Power Sector, User Guide Version 20.0
# (PDF PAGE 4): "Table S: Weighted average emission factor, simple
# operating margin (OM), build margin (BM) and combined margin (CM) of the Indian Grid for FY 2023-24 (adjusted for
# cross-border electricity transfers) (including RES & Captive power injection into grid), in tCO2/MWh" /
# "Average OM BM CM" / "0.727  0.962  0.552  0.757".
# Transmission and distribution losses are NOT added.                                    [ASSUMPTION]
GRID_INDIA       = 0.727  # t CO2/MWh (= 727 g CO2/kWh)                                           [VERIFIED]
# 0.710 = FY 2024-25. When this script was run, this value had not been checked in a CEA document, so the
# label written to route_sensitivity.csv says UNVERIFIED. CEA User Guide Version 22.0 (Appendix C, Table B)
# later confirmed 0.710 for FY 2024-25 and 0.727 for FY 2023-24. Used only as a sensitivity case.
GRID_INDIA_V21_REPORTED = 0.710
GRID_RENEWABLE   = 0.0    # t CO2/MWh, operational emissions only (life-cycle of PV/wind excluded) [ASSUMPTION]

# Scrap shares evaluated: study scenarios                                                 [ASSUMPTION]
SCRAP_SHARES     = [0.0, 0.50, 0.95]
GRIDS = [("Renewable", GRID_RENEWABLE), ("India grid", GRID_INDIA)]

# Plant data: CCTS iron & steel GHG Emission Intensity schedule, baseline year 2023-24,
# GEI in tCO2e per tonne of EQUIVALENT PRODUCT, output in tonnes of equivalent product.
# MAIN: "Greenhouse Gases Emission Intensity Target (Amendment) Rules, 2025" - DRAFT, MoEFCC, G.S.R. 517(E),
#   26 June 2026 (MoEFCC_2026_517E_pdf.txt). English schedule, "Sector: Iron and Steel", column heads
#   "Registration Number", "(Tonne)", "(tCO2 equivalent /Tonne of equivalent product)". Rows read:
#   PDF PAGE 31: "Bhilai Steel Plant, Bhilai, Durg, Chhattisgarh, 490001" / "INSOE0005CG 5675148 3.1487 - 2.9551"
#   PDF PAGE 32: "Tata Steel Ltd, Jamshedpur, ..." / "INSOE0006JH 10703390 2.3630 - 2.2540"
#   PDF PAGE 32: "JSW STEEL LTD, Vijayanagar Works, ..." / "INSOE0008KA 11739649 2.6178 - 2.4840"
#   PDF PAGE 32: "JSW Steel LTD, Geetapuram, Dolvi, ..." / "INSOE0011MH 8904781 2.7222 - 2.5775"   [VERIFIED]
# SENSITIVITY: original draft, G.S.R. 405(E), 23 June 2025 (MoEFCC_2025_405E_pdf.txt), PDF PAGE 72:
#   "5.   Bhilai Steel Plant, Bhilai, Durg," / "... 5675148 3.1487 3.0713 2.9553"
#   "6.   Tata Steel Ltd, Jamshedpur," / "... 10703390 2.3804 2.3362 2.2699"
#   "8.   JSW STEEL LTD, Vijayanagar" / "... 11739649 2.6293 2.5754 2.4945"
#   "Dolvi, Raigad, Maharastra," / "... 8904781 2.6662 2.6107 2.5275"                              [VERIFIED]
# Plant choice [ASSUMPTION]: the four ore-based plants that have a rebar EPD in Step II.
PLANTS = [
    # name,                  CCTS reg. no.,  output_2026,  GEI_2026,  output_2025,  GEI_2025
    ("Tata Steel Jamshedpur", "INSOE0006JH", 10703390, 2.3630, 10703390, 2.3804),
    ("JSW Steel Vijayanagar", "INSOE0008KA", 11739649, 2.6178, 11739649, 2.6293),
    ("JSW Steel Dolvi",       "INSOE0011MH",  8904781, 2.7222,  8904781, 2.6662),
    ("SAIL Bhilai",           "INSOE0005CG",  5675148, 3.1487,  5675148, 3.1487),
]
SHORT = ["Jamshedpur", "Vijayanagar", "Dolvi", "Bhilai"]
for p, key in zip(PLANTS, SHORT):
    assert key in p[0]
    assert 1.0 < p[3] < 5.0 and 1.0 < p[5] < 5.0, f"GEI out of range (unit error?): {p[0]}"
    assert p[2] > 1e6, f"output out of range (unit error?): {p[0]}"

with open(f"{OUT}/plant_inputs.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["plant", "ccts_registration", "output_2023_24_t_rev2026", "GEI_2023_24_rev2026",
                 "output_2023_24_t_draft2025", "GEI_2023_24_draft2025"])
    for p in PLANTS:
        wr.writerow([p[0], p[1], p[2], f"{p[3]:.4f}", p[4], f"{p[5]:.4f}"])

# Step II evidence, read ONLY for the side-by-side table of PART 5. It enters no calculation.
EST_CSV = "outputs_phase2/route_estimates.csv"
INP_CSV = "outputs_phase2/epd_inputs.csv"
for _need in (EST_CSV, INP_CSV):
    if not os.path.exists(_need):
        raise FileNotFoundError("Run phase2_route_estimates_colab.py first: " + _need + " not found.")
GEI = [p[3] for p in PLANTS]                # t CO2e / t equivalent product, 2026 revised draft (main)
GEI_2025 = [p[5] for p in PLANTS]           # 2025 draft (sensitivity)

# =============================================================================
# PART 2 -- MODEL
# =============================================================================
def sec(s):                 # MWh/t; linear between Vogl's two end points          [ASSUMPTION]
    return SEC_H2DRI_0SCRAP - s * (SEC_H2DRI_0SCRAP - SEC_EAF_100SCRAP)

def intensity(s, g):        # t CO2 / t liquid steel
    return DIRECT_CO2 + sec(s) * g + REST_FUEL_GJ * EF_NATGAS + REST_ELEC_MWH * g

def route(s, g, fuel_gj=REST_FUEL_GJ, rest_mwh=REST_ELEC_MWH):   # same model with the two "rest" terms as arguments
    return DIRECT_CO2 + sec(s) * g + fuel_gj * EF_NATGAS + rest_mwh * g

def breakeven_grid(I_target, s, fuel_gj=REST_FUEL_GJ, rest_mwh=REST_ELEC_MWH):
    """g CO2/kWh at which route(s, g) = I_target."""
    return (I_target - DIRECT_CO2 - fuel_gj * EF_NATGAS) / (sec(s) + rest_mwh) * 1000

def scrap_at_intensity(I_target, g, fuel_gj=REST_FUEL_GJ, rest_mwh=REST_ELEC_MWH):
    """Scrap share s where route(s, g) = I_target (analytic inverse; can fall outside 0..1)."""
    return (SEC_H2DRI_0SCRAP + rest_mwh - (I_target - DIRECT_CO2 - fuel_gj * EF_NATGAS) / g) \
           / (SEC_H2DRI_0SCRAP - SEC_EAF_100SCRAP)

def share_text(sx):
    """Break-even share as text: a number inside 0..1, else which side it falls on."""
    if sx < 0:
        return "below 0 (route is lower than the plant at every scrap share)"
    if sx > 1:
        return "above 1 (route is higher than the plant at every scrap share)"
    return z0(sx)

# Regression guards (the only two): they stop the run if a route constant is edited by mistake.
assert round(intensity(0.0, 0.727), 6) == 2.808160, "regression guard: I(0, 0.727) changed; were the route constants edited?"
assert round(intensity(0.95, 0.727), 6) == 0.865362, "regression guard: I(0.95, 0.727) changed; were the route constants edited?"
# Property checks
assert abs(sec(0.0) - SEC_H2DRI_0SCRAP) < 1e-12 and abs(sec(1.0) - SEC_EAF_100SCRAP) < 1e-12
assert all(abs(route(s, g) - intensity(s, g)) < 1e-15 for s in SCRAP_SHARES for _, g in GRIDS)
assert intensity(0.0, GRID_INDIA) > intensity(0.5, GRID_INDIA) > intensity(0.95, GRID_INDIA)   # more scrap -> lower on the grid
assert len({round(intensity(s, 0.0), 12) for s in SCRAP_SHARES}) == 1                          # renewable: no effect of scrap

# =============================================================================
# PART 3 -- SCENARIO TABLE
# =============================================================================
print("=" * 78)
print(f"{'scrap':>6} {'electricity':>12} {'g/kWh':>7} {'SEC MWh/t':>10} {'I_new t/t':>10}")
print("=" * 78)
scen = []
for s in SCRAP_SHARES:
    for name, g in GRIDS:
        scen.append((s, name, g, sec(s), intensity(s, g)))
        print(f"{s:6.2f} {name:>12} {g*1000:7.1f} {sec(s):10.3f} {intensity(s, g):10.3f}")
with open(f"{OUT}/scenario_table.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["scrap_share", "electricity", "grid_gCO2_per_kWh", "SEC_MWh_per_t", "I_new_tCO2_per_t"])
    for s, name, g, se, I in scen:
        wr.writerow([s, name, f"{g*1000:.6f}", f"{se:.6f}", f"{I:.6f}"])

# =============================================================================
# PART 4 -- BREAK-EVEN PER PLANT (own baseline)
# =============================================================================
print("\nBreak-even per plant (route = plant's own GEI):")
print(f"  {'plant':24s} {'GEI':>7} {'scrap share on India grid':>28}   grid factor g/kWh at s=0 / 0.50 / 0.95")
be_rows = []
for version, gei_list in (("rev2026 (main)", GEI), ("draft2025 (sensitivity)", GEI_2025)):
    for p, gei in zip(PLANTS, gei_list):
        sx = scrap_at_intensity(gei, GRID_INDIA)
        if 0.0 <= sx <= 1.0:      # inside the range the break-even must reproduce the plant's GEI
            assert abs(intensity(sx, GRID_INDIA) - gei) < 1e-12
        else:                     # outside the range the route never crosses the plant's GEI for 0 <= s <= 1
            lo, hi = intensity(1.0, GRID_INDIA), intensity(0.0, GRID_INDIA)
            assert gei > hi if sx < 0 else gei < lo
        txt = share_text(sx)
        bg = [breakeven_grid(gei, s) for s in SCRAP_SHARES]
        for s, b in zip(SCRAP_SHARES, bg):
            assert abs(intensity(s, b / 1000) - gei) < 1e-12
        be_rows.append([p[0], version, f"{gei:.4f}", z0(sx), txt] + [f"{b:.3f}" for b in bg])
        if version.startswith("rev2026"):
            print(f"  {p[0]:24s} {gei:7.4f} {(f'{100*sx:.1f}%' if 0 <= sx <= 1 else txt.split(' (')[0]):>28}   "
                  + " / ".join(f"{b:7.1f}" for b in bg))
with open(f"{OUT}/plant_breakeven.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["plant", "GEI_version", "GEI_tCO2e_per_t", "breakeven_scrap_share_analytic_india_grid",
                 "breakeven_scrap_share_india_grid_text"]
                + [f"breakeven_grid_gCO2_per_kWh_scrap{int(round(s*100))}" for s in SCRAP_SHARES])
    wr.writerows(be_rows)

# =============================================================================
# PART 5 -- ROUTE MODEL AND EPD EVIDENCE, SIDE BY SIDE (NOT COMBINED)
#   The two blocks of this table are DIFFERENT QUANTITIES:
#     route model : t CO2 per t liquid steel; on-site emissions and electricity only (boundary of Vogl et al.)
#     EPD values  : kg CO2e per kg rebar; cradle to gate (modules A1-A3), declared by the producers
#   They are shown together ONLY to see whether the direction of the model with scrap share agrees with
#   declared values. No ratio is taken and nothing is scaled. Do not subtract or divide one block by the other.
# =============================================================================
NOTE_MODEL = ("route model; t CO2 per t steel; on-site and electricity only; a different quantity from the EPD rows; "
              "shown only to compare the direction with scrap share; no ratio, no scaling")
NOTE_EPD = ("declared EPD values; kg CO2e per kg rebar; cradle to gate (A1-A3); a different quantity from the route "
            "model rows; shown only to compare the direction with scrap share; no ratio, no scaling")
with open(EST_CSV, encoding="utf-8") as f:
    _est = {r["route"]: r for r in csv.DictReader(f)}
with open(INP_CSV, encoding="utf-8") as f:
    _inp = [r for r in csv.DictReader(f) if not r["input"].startswith("PRIOR")]
ev_rows = []
for name, g in GRIDS:
    for s in SCRAP_SHARES:
        ev_rows.append(["route model", f"{name}, scrap share {s:.2f}", name, f"{s:.2f}", "", "route intensity",
                        f"{intensity(s, g):.6f}", "t CO2 per t liquid steel", NOTE_MODEL])
for grp in ("ore-based", "scrap"):
    if grp not in _est:
        raise KeyError(f"{EST_CSV} has no row for the route group '{grp}'")
    vals = [float(r["gwp_total_a1a3_kg_per_kg"]) for r in _inp if r["route_group"] == grp]
    if not vals:
        raise KeyError(f"{INP_CSV} has no EPD row for the route group '{grp}'")
    e = _est[grp]
    # the range in the two Step II files must be the same (two files, one set of inputs)
    assert abs(min(vals) - float(e["input_min"])) < 5e-7 and abs(max(vals) - float(e["input_max"])) < 5e-7, \
        f"input range of the {grp} group differs between {INP_CSV} and {EST_CSV}"
    assert len(vals) == int(e["n"]), f"number of EPD inputs of the {grp} group differs between the two Step II files"
    rej = e["common_mean_rejected_at_5pct"]
    unit = "kg CO2e per kg rebar"
    for item, v in (("smallest EPD input", min(vals)), ("largest EPD input", max(vals)),
                    ("estimate, fixed effect with prior" + (" (common mean rejected: do not quote its interval)"
                                                             if rej == "yes" else ""), float(e["mu_post"])),
                    ("estimate, random effects, EPDs only", float(e["mu_random_effects"]))):
        ev_rows.append(["Step II EPD evidence", f"{grp} group ({len(vals)} EPD values)", "", "", grp, item,
                        f"{v:.6f}", unit, NOTE_EPD])
with open(f"{OUT}/route_vs_epd_evidence.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["block", "row", "electricity", "scrap_share", "epd_route_group", "item", "value", "unit", "note"])
    wr.writerows(ev_rows)
print("\nRoute model and EPD evidence, side by side (different quantities; no ratio, no scaling):")
for r in ev_rows:
    print(f"  {r[0]:21s} {r[1]:34s} {r[5][:46]:46s} {r[6]:>9}  {r[7]}")
_model_falls = intensity(SCRAP_SHARES[0], GRID_INDIA) > intensity(SCRAP_SHARES[-1], GRID_INDIA)
_epd_falls = float(_est["scrap"]["input_max"]) < float(_est["ore-based"]["input_min"])
print(f"  Direction only: model on the India grid falls with scrap share: {_model_falls}; every scrap-group EPD value "
      f"is below every ore-based EPD value: {_epd_falls}. With renewable electricity the model does not change with "
      f"scrap share.")

# =============================================================================
# PART 6 -- ROUTE SENSITIVITY: grid factor and system boundary
# =============================================================================
SENS_CASES = [
    ("main (CEA Version 20.0, FY 2023-24; with pelletizing and rolling)", GRID_INDIA, REST_FUEL_GJ, REST_ELEC_MWH),
    ("grid 0.710 (reported for CEA Version 21.0, FY 2024-25; document not read; UNVERIFIED, sensitivity only)",
     GRID_INDIA_V21_REPORTED, REST_FUEL_GJ, REST_ELEC_MWH),
    ("main grid, without pelletizing and rolling terms (boundary of the Vogl et al. energy figure)", GRID_INDIA, 0.0, 0.0),
]
print("\nRoute sensitivity (grid factor and boundary); break-even scrap share per plant:")
with open(f"{OUT}/route_sensitivity.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["case", "grid_tCO2_per_MWh"] + [f"I_scrap{int(round(s*100))}_tCO2_per_t" for s in SCRAP_SHARES]
                + [f"breakeven_scrap_share_{k}" for k in SHORT])
    for name, g, fu, re_ in SENS_CASES:
        Is = [route(s, g, fu, re_) for s in SCRAP_SHARES]
        sx = [scrap_at_intensity(gei, g, fu, re_) for gei in GEI]
        for x, gei in zip(sx, GEI):
            if 0.0 <= x <= 1.0:
                assert abs(route(x, g, fu, re_) - gei) < 1e-12
        wr.writerow([name, f"{g:.3f}"] + [f"{v:.6f}" for v in Is] + [share_text(x) for x in sx])
        print(f"  {name[:58]:58s} I = " + " / ".join(f"{v:.3f}" for v in Is) + "   s* = "
              + ", ".join((f"{100*x:.1f}%" if 0 <= x <= 1 else share_text(x).split(' (')[0]) for x in sx))

# =============================================================================
# PART 7 -- SWEEPS (CSV)
# =============================================================================
s_grid = np.linspace(0, 1, 101)
g_grid = np.linspace(0, 0.800, 161)          # 0-800 g/kWh
with open(f"{OUT}/sweep_scrap.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["scrap_share", "I_renewable", "I_india_grid"])
    for s in s_grid:
        wr.writerow([f"{s:.2f}", f"{intensity(s, 0):.6f}", f"{intensity(s, GRID_INDIA):.6f}"])
with open(f"{OUT}/sweep_grid.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f); wr.writerow(["grid_gCO2_per_kWh"] + [f"I_s{int(round(s*100))}" for s in SCRAP_SHARES])
    for g in g_grid:
        wr.writerow([f"{g*1000:.1f}"] + [f"{intensity(s, g):.6f}" for s in SCRAP_SHARES])

# =============================================================================
# PART 8 -- FIGURES (1000 dpi + PDF TIFF, LZW)
# =============================================================================
def save(fig, path, dpi=1000):
    fig.tight_layout()
    fig.savefig(path, format="tiff", dpi=dpi, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(path.rsplit(".", 1)[0] + ".pdf", format="pdf", bbox_inches="tight")   # vector copy
    plt.close(fig); print("Wrote", path)

def clean(ax):
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

COLS = ["#2c7fb8", "#807dba", "#f16913"]                       # scrap 0, 50, 95%
PCOL = ["#111111", "#666666", "#a6611a", "#c51b7d"]            # plants J, V, D, B
PLS = ["--", ":", "-.", (0, (5, 1, 1, 1, 1, 1))]

# Fig 1 -- scrap-share sweep with the four plant GEI lines
fig, ax = plt.subplots(figsize=(8, 4.8), dpi=1000)
ax.plot(s_grid * 100, [intensity(s, GRID_INDIA) for s in s_grid], color="#f16913", lw=2,
        label=f"H$_2$-DRI/EAF, India grid ({GRID_INDIA*1000:.0f} g/kWh)")
ax.plot(s_grid * 100, [intensity(s, 0) for s in s_grid], color="#41ab5d", lw=2, label="H$_2$-DRI/EAF, renewable (operational 0 g/kWh)")
for k, gei, c, ls in zip(SHORT, GEI, PCOL, PLS):
    ax.axhline(gei, color=c, ls=ls, lw=1.2, label=f"{k} GEI {gei:.4f}")
    sx = scrap_at_intensity(gei, GRID_INDIA)
    if 0 <= sx <= 1:
        ax.plot(sx * 100, gei, "D", color=c, ms=5, zorder=5)
        ax.annotate(f"{100*sx:.1f}%", (sx * 100, gei), textcoords="offset points", xytext=(6, 5), fontsize=8, color=c,
                    bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none"))
ax.set_xlabel("Scrap share in EAF charge (%)", fontsize=12)
ax.set_ylabel(r"t CO$_2$ / t steel", fontsize=12)
ax.set_title("Route intensity vs scrap share and plant baselines", fontsize=14, fontweight="bold", pad=15)
ax.set_xlim(0, 100); ax.set_ylim(0, None); clean(ax)
ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=3)
save(fig, f"{OUT}/fig_scrap_sweep.tiff")

# Fig 2 -- grid-intensity sweep with the four plant GEI lines
fig, ax = plt.subplots(figsize=(8, 4.8), dpi=1000)
for s, c in zip(SCRAP_SHARES, COLS):
    ax.plot(g_grid * 1000, [intensity(s, g) for g in g_grid], color=c, lw=2, label=f"H$_2$-DRI/EAF, scrap {int(round(s*100))}%")
for k, gei, c, ls in zip(SHORT, GEI, PCOL, PLS):
    ax.axhline(gei, color=c, ls=ls, lw=1.2, label=f"{k} GEI {gei:.4f}")
ax.axvline(GRID_INDIA * 1000, color="#999999", ls="-", lw=1.0, label=f"India grid {GRID_INDIA*1000:.0f} g/kWh (CEA, FY 2023-24)")
for s, c in zip(SCRAP_SHARES, COLS):          # break-even markers where a route line crosses a plant line
    for gei in GEI:
        be = breakeven_grid(gei, s)
        if be <= 800:
            ax.plot(be, gei, "o", color=c, ms=4, mec="black", mew=0.5, zorder=5)
ax.set_xlabel(r"Grid emission factor (g CO$_2$ / kWh)", fontsize=12)
ax.set_ylabel(r"t CO$_2$ / t steel", fontsize=12)
ax.set_title("Route intensity vs electricity carbon intensity", fontsize=14, fontweight="bold", pad=15)
ax.set_xlim(0, 800); ax.set_ylim(0, None); clean(ax)
ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=3)
save(fig, f"{OUT}/fig_grid_sweep.tiff")

import platform
with open(f"{OUT}/software_versions_phase3_route_model_colab.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["package", "version"])
    wr.writerows([["python", platform.python_version()], ["numpy", np.__version__],
                  ["matplotlib", matplotlib.__version__]])

print("\nStep III done. Files in", OUT)
