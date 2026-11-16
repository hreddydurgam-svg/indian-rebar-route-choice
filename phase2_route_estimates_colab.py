# -*- coding: utf-8 -*-
"""
phase2_route_estimates_colab.py
===============================
Step II - Two route estimates for Indian steel rebar from rebar-only Environmental Product
Declarations (EPDs): 1 kg rebar, cradle-to-gate (modules A1-A3), total climate indicator
(GWP-total, EN 15804+A2), kg CO2e / kg rebar.

  Route group "ore-based": rebar from ore-based steel plants (4 values).
  Route group "scrap"    : rebar from scrap melted in induction or electric arc furnaces (3 values).
Each group gets its own normal-normal update with the same global prior, and is also reported
without the prior. A spread test over all seven values shows with a number that one common value
for Indian rebar does not fit. The same test is run INSIDE each group, and a random-effects
estimate (DerSimonian-Laird) is given beside the fixed-effect one.

HOW TO QUOTE THE RESULT
  * Where the test rejects one common mean inside a group at the 5% level (column
    common_mean_rejected_at_5pct = yes; in this data set: the scrap group), the fixed-effect 90%
    interval (ci90_lower, ci90_upper) MUST NOT be quoted: it is the interval of a common mean that the
    data reject. Quote the input range (input_min to input_max) and the random-effects interval instead.
  * plant_epd_values.csv lists each plant's EPD value for the table of the paper that shows it beside
    the plant's gazette intensity (GEI). The two are NEVER combined (no ratio, no scaling) in this or any
    later script: an EPD (cradle to gate, product level, CO2e) and a GEI (plant level, per tonne of
    equivalent product) come from different accounting systems.

HOW TO RUN IN GOOGLE COLAB
  1. Upload this file (Files panel, left side) or paste it into ONE code cell.
  2. Run:   !python phase2_route_estimates_colab.py      (or run the cell)
  3. Then run phase3_route_model_colab.py in the SAME runtime (it reads route_estimates.csv and epd_inputs.csv).
Outputs are written to ./outputs_phase2/
  epd_inputs.csv               the seven values used, with assumed SD, weight, route as the EPD states it
                               (route_note), source file and page
  route_estimates.csv          route, n, mu_post, sd_post, ci90_lower, ci90_upper (main case, fixed effect, with
                               prior); then input_min, input_max, Cochran Q, df, p-value, I^2,
                               common_mean_rejected_at_5pct, and the random-effects estimate with its 90% interval
  route_sensitivity_cases.csv  per route: prior SD 10/30/50%, EPD SD 5/20%, equal EPD SD, leave one out, no prior
  heterogeneity.csv            Cochran Q and I^2 over all seven values and within each group
  plant_epd_values.csv         plant, epd_kgCO2e_per_kg (the four ore-based plants; for the paper's table only)
  epd_screening.csv            every Indian rebar-related EPD found, with the decision and the reason
  audit_route_unrounded.csv    every weight and contribution, unrounded
  software_versions_phase2_route_estimates_colab.csv   python, scipy, matplotlib versions of the run
  fig_route_estimates.tiff     1000 dpi + PDF, LZW

WHAT THIS SCRIPT DOES NOT SHOW (state these limits in the paper)
  * Seven values from six EPD documents of six producers (Tata Steel, JSW Steel, SAIL, ARS Steels,
    Prakash Ferrous, Jindal Steel & Power). The two JSW values come from ONE document (two sites).
    The producers chose to publish these EPDs; this is not a random sample of Indian rebar.
  * EPDs are self-declared by the producers. They differ
    in data year (April 2019 - March 2020 up to FY 2024-25), background database and scrap rules.
  * The standard deviations (30% for the prior, 10% for each EPD) are ASSUMPTIONS. No EPD prints an
    uncertainty. The sensitivity table shows what changes when they change.
  * The fixed-effect 90% interval is the interval for ONE COMMON MEAN within a route under those assumptions.
    It is NOT the plant-to-plant spread, which is wider (see heterogeneity.csv and the input values).
  * The random-effects interval uses a normal quantile. With three or four values the between-EPD variance
    is poorly estimated, so this interval is also likely too narrow. No small-sample correction is applied.
  * With an SD proportional to the value, low values get more weight. The "equal EPD SD" case removes this.
  * The scrap group has only three values (two without the special case). One of the three (Jindal) is a
    special case: electric arc furnace with 90% renewable electricity. The scrap estimate is weak evidence.
  * Rebar made from coal-based sponge iron in induction furnaces, a large part of Indian rebar supply, has
    no usable EPD in this set. Therefore NO national value for Indian rebar is given by this script.
  * The JSW Dolvi value is placed in the ore-based group (ore-based plant with blast furnace and sponge
    iron plant), but its EPD states the CONARC process (hot metal and sponge iron in an electric arc furnace),
    not a basic oxygen furnace.
  * The Tata Steel EPD covers rebar from Jamshedpur steel rolled at Jamshedpur mills AND at subsidiaries and
    contract mills. It is used as the value for the Jamshedpur plant.                            [ASSUMPTION]

Every input below carries a tag:
  [VERIFIED]   value read in the named source text (page marker and quote given)
  [ASSUMPTION] study choice; must be stated as an assumption in the paper
"""

# =============================================================================
# PART 0 -- COLAB DEPENDENCIES (runs only inside Colab; skipped elsewhere)
# =============================================================================
import sys, subprocess, glob, os

IN_COLAB = "google.colab" in sys.modules or os.path.isdir("/content")
if IN_COLAB:
    cmds = [
        # accept the Microsoft core-fonts EULA non-interactively, then install Arial
        "echo ttf-mscorefonts-installer msttcorefonts/accept-mscorefonts-eula select true | sudo debconf-set-selections",
        "sudo apt-get -qq update -y",
        "sudo DEBIAN_FRONTEND=noninteractive apt-get -qq install -y ttf-mscorefonts-installer",
        "rm -rf ~/.cache/matplotlib",
        f"{sys.executable} -m pip install -q matplotlib pillow numpy scipy",
    ]
    for c in cmds:
        subprocess.run(c, shell=True, check=False)

# =============================================================================
# IMPORTS
# =============================================================================
from decimal import Decimal, getcontext
import csv
from scipy.stats import chi2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

getcontext().prec = 50

# register Arial if installed; otherwise fall back silently to DejaVu Sans
for f in glob.glob("/usr/share/fonts/truetype/msttcorefonts/*.ttf"):
    try:
        fm.fontManager.addfont(f)
    except Exception:
        pass
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]

OUT = "outputs_phase2"
os.makedirs(OUT, exist_ok=True)
SEED = 12345          # not used: this script has no random step (kept for the common script layout)

def dec_str(d): return format(d, "f")
def f6(d): return str(Decimal(d).quantize(Decimal("0.000001")))
def r3(d): return Decimal(d).quantize(Decimal("0.001"))
def r1(d): return Decimal(d).quantize(Decimal("0.1"))

# exact two-sided 90% normal quantile
Z90 = Decimal("1.6448536269514722")

# =============================================================================
# PART 1 -- THE ESTIMATOR AND ITS TEST
#   mu_post = (sum w_i x_i) / (sum w_i),  w_i = 1/sigma_i^2,  sigma_post = 1/sqrt(sum w_i)
# =============================================================================
def pooled(mu0, sd0, obs):
    """Normal-normal update. obs = list of (name, value, sd). mu0 = None -> no prior (EPDs only)."""
    names = [n for n, _, _ in obs]
    x = [v for _, v, _ in obs]
    s = [sd for _, _, sd in obs]
    if mu0 is not None:
        names, x, s = ["PRIOR"] + names, [mu0] + x, [sd0] + s
    w = [Decimal(1) / (si * si) for si in s]
    W = sum(w)
    S = sum(wi * xi for wi, xi in zip(w, x))
    mu = S / W
    sd = Decimal(1) / W.sqrt()
    return {"names": names, "x": x, "s": s, "w": w, "W": W, "S": S, "mu": mu, "sd": sd,
            "lo": mu - Z90 * sd, "hi": mu + Z90 * sd, "weights": [wi / W for wi in w]}

# test 1: with one observation of equal SD the mean is the midpoint, and the weights sum to 1
_r = pooled(Decimal("1"), Decimal("0.1"), [("a", Decimal("3"), Decimal("0.1"))])
assert abs(_r["mu"] - Decimal("2")) < Decimal("1e-40")
assert abs(sum(_r["weights"]) - 1) < Decimal("1e-40")
# test 2: without a prior and with equal SDs the result is the plain average, SD = sd / sqrt(n)
_r = pooled(None, None, [("a", Decimal("1"), Decimal("0.2")), ("b", Decimal("2"), Decimal("0.2")),
                            ("c", Decimal("6"), Decimal("0.2")), ("d", Decimal("7"), Decimal("0.2"))])
assert abs(_r["mu"] - Decimal("4")) < Decimal("1e-40") and abs(_r["sd"] - Decimal("0.1")) < Decimal("1e-40")

def cochran(obs):
    """Cochran Q = sum w_i (x_i - mean_w)^2, w_i = 1/sd_i^2; df = k - 1; I^2 = max(0, (Q - df)/Q)."""
    w = [Decimal(1) / (sd * sd) for _, _, sd in obs]
    x = [v for _, v, _ in obs]
    mean_w = sum(wi * xi for wi, xi in zip(w, x)) / sum(w)
    Q = sum(wi * (xi - mean_w) ** 2 for wi, xi in zip(w, x))
    df = len(obs) - 1
    I2 = max(Decimal(0), (Q - df) / Q) if Q > 0 else Decimal(0)
    return {"k": len(obs), "mean_w": mean_w, "Q": Q, "df": df, "I2": I2, "p": float(chi2.sf(float(Q), df)),
            "min": min(x), "max": max(x)}

# test 3: identical values give Q = 0 and I^2 = 0
_c = cochran([("a", Decimal("2"), Decimal("0.2")), ("b", Decimal("2"), Decimal("0.3"))])
assert _c["Q"] == 0 and _c["I2"] == 0

def random_effects(obs):
    """Random-effects estimate of a group, DerSimonian-Laird method, EPD inputs only (NO prior).
         w_i   = 1 / sd_i^2
         tau^2 = max(0, (Q - df) / (sum w - sum w^2 / sum w))      between-EPD variance
         w*_i  = 1 / (sd_i^2 + tau^2);  mu = sum w*_i x_i / sum w*_i;  sd = 1 / sqrt(sum w*_i)
       90% interval = mu +/- 1.645 sd (normal quantile, no small-sample correction).
       [ASSUMPTION] choice of the DerSimonian-Laird method; other tau^2 estimators give other intervals."""
    x = [v for _, v, _ in obs]
    v = [sd * sd for _, _, sd in obs]
    w = [Decimal(1) / vi for vi in v]
    c = cochran(obs)
    denom = sum(w) - sum(wi * wi for wi in w) / sum(w)
    tau2 = max(Decimal(0), (c["Q"] - c["df"]) / denom) if denom > 0 else Decimal(0)
    ws = [Decimal(1) / (vi + tau2) for vi in v]
    mu = sum(wi * xi for wi, xi in zip(ws, x)) / sum(ws)
    sd = Decimal(1) / sum(ws).sqrt()
    return {"tau2": tau2, "mu": mu, "sd": sd, "lo": mu - Z90 * sd, "hi": mu + Z90 * sd}

# test 4: when Q <= df the between-EPD variance is zero and the result equals the fixed-effect one (no prior)
_o = [("a", Decimal("2.0"), Decimal("0.5")), ("b", Decimal("2.1"), Decimal("0.5")), ("c", Decimal("1.9"), Decimal("0.4"))]
_re, _fe = random_effects(_o), pooled(None, None, _o)
assert _re["tau2"] == 0 and abs(_re["mu"] - _fe["mu"]) < Decimal("1e-40") and abs(_re["sd"] - _fe["sd"]) < Decimal("1e-40")
# test 5: two far-apart values with equal SD: tau^2 > 0, mean = midpoint, SD larger than the fixed-effect SD
_o = [("a", Decimal("1"), Decimal("0.1")), ("b", Decimal("3"), Decimal("0.1"))]
_re, _fe = random_effects(_o), pooled(None, None, _o)
assert _re["tau2"] > 0 and abs(_re["mu"] - Decimal("2")) < Decimal("1e-40") and _re["sd"] > _fe["sd"]

# =============================================================================
# PART 2 -- INPUTS (single source of truth)
# =============================================================================
# Prior: worldsteel LCI data for steel products, 2022 data release (provided 18-08-2025),
# 1 kg global rebar, cradle-to-gate, excl. end-of-life recycling credit,
# "IPCC AR6 GWP 100, incl biogenic CO2" = 1.73 kg CO2e/kg
#                                                 [VERIFIED by the authors against the worldsteel data report]
# (the worldsteel report is not shared; only this aggregated value is used, under worldsteel terms)
# The same prior is used for both route groups. It is a global average over routes, so it is not
# route-specific; each route is therefore also reported without it.                       [ASSUMPTION]
MU0 = Decimal("1.73")
REL_SD0 = Decimal("0.30")                                                 # [ASSUMPTION]
SIGMA0 = REL_SD0 * MU0
REL_SD_EPD = Decimal("0.10")                                              # [ASSUMPTION] no EPD prints an uncertainty

# Evidence: Indian rebar EPDs, modules A1-A3, GWP-total, printed per tonne (kg CO2e/t), entered here as
# printed value / 1000 = kg CO2e/kg. Source texts carry "===== PDF PAGE n =====" markers.
# Fields: key, route group, company, product, registration, value, special case?, plant (or None),
#         route as stated, data period, source file, page of the value.
EPD = [
    # --- Tata Steel, EPD S-P-05018 (PDF PAGE 1: "EPD registration number:  S-P-05018") -------------------
    # PDF PAGE 15, table row " Climate Change kg CO2 eq. 2.52E+03 2.16E+00 2.05E+02 6.10E+00 ..." under the
    #   header " Parameter Unit A1 A2 A3 C1 C2 C3 C4 D". The table has no A1-A3 column, so
    #   A1 + A2 + A3 = 2520 + 2.16 + 205 = 2727.16 kg/t.                  [VERIFIED]
    #   PDF PAGE 17 prints the rounded total: "The total cradle to gate impact is 2.73E+03 kg CO2 eq."
    # Route, PDF PAGE 9: "represents 100% Blast Furnace (BF) with Basic Oxygen Furnace (BOF) route".
    # Period, PDF PAGE 5: "The final results mentioned in this study, represents the production tonnage
    #   manufactured during the year April 2020 to March 2021" (primary data FY 2018 to FY 2021).
    # Scope, PDF PAGE 5: "manufactured by Tata Steel Limited at its Jamshedpur facility as well as from its
    #   subsidiaries and external steel processing centres in India".
    # The brand name "Tiscon" is in the file name only: the cover page text is not machine-readable.
    ("Tata_Steel_rebar_S-P-05018", "ore-based", "Tata Steel Limited", "Steel rebar", "S-P-05018",
     (Decimal("2520") + Decimal("2.16") + Decimal("205")) / 1000, False, "Tata Steel Jamshedpur",
     "100% blast furnace with basic oxygen furnace (BF-BOF)", "April 2020 - March 2021 (primary data FY 2018 to FY 2021)",
     "Tata_Tiscon_Environmental_Product_Declaration_Final_pdf.txt", "15"),
    # --- JSW Steel, EPD S-P-04326 (PDF PAGE 1: "S-P-04326 2021-10-21 2026-10-20 India"); ONE document, two sites
    # PDF PAGE 15, "Vijayanagar", "Table 11  Environmental impacts for 1- tonne TMT Bar",
    #   "Parameter Unit A1-A3 C1 C2 C3 C4 D" / "GWP - total kg CO2 eq. 2.42E+03 0.00E+00 ..."      [VERIFIED]
    # Route, PDF PAGE 7: "Products are manufactured using the Blast Furnace with Basic Oxygen Furnace (BF+BOF)
    #   route"; PDF PAGE 13: "Vijayanagar- Steel is produced predominantly by two process routes: the blast
    #   furnace/basic oxygen furnace route and the electric arc furnace route ... COREX technology is also used".
    # Period, PDF PAGE 4: "The reference period for the primary data used within this EPD is the April 2019 to
    #   March 20 20."
    ("JSW_Vijayanagar_TMT_S-P-04326", "ore-based", "JSW Steel Limited", "TMT bar, Vijayanagar", "S-P-04326",
     Decimal("2420") / 1000, False, "JSW Steel Vijayanagar",
     "blast furnace and COREX with basic oxygen furnace; electric arc furnace route also on site",
     "April 2019 - March 2020", "Final_EPD_Long_Finished_Steel_products_JSW_Steel_pdf.txt", "15"),
    # PDF PAGE 18, "Dolvi", "Table 19 Environmental impacts for 1- tonne TMT Bar",
    #   "GWP - total kg CO2 eq. 3.25E+03 0.00E+00 4.13E-03 ..."                                    [VERIFIED]
    # Route, PDF PAGE 13: "Dolvi - Steel is produced by CONARC process."; PDF PAGE 8: "The Dolvi unit features
    #   a sponge iron plant, a blast furnace, pellet plant, coke-oven plant, sinter plants" and "The hot iron
    #   produced in the Blast Furnace and Direct Reduced Iron (DRI) is transferred into SMS (Steel Melting Shop)
    #   and CONARC to produce liquid steel".
    # Group choice: plant with its own ore-based ironmaking, so "ore-based"; it is NOT a BF-BOF plant. [ASSUMPTION]
    ("JSW_Dolvi_TMT_S-P-04326", "ore-based", "JSW Steel Limited", "TMT bar, Dolvi", "S-P-04326",
     Decimal("3250") / 1000, False, "JSW Steel Dolvi",
     "CONARC process (hot metal and sponge iron in an electric furnace)",
     "April 2019 - March 2020", "Final_EPD_Long_Finished_Steel_products_JSW_Steel_pdf.txt", "18"),
    # --- SAIL Bhilai Steel Plant, EPD-IES-0031470 (PDF PAGE 1: "EPD registration number: EPD-IES-0031470") -
    # PDF PAGE 14, "Indicator Unit A1-A3 A4 C1 C2 C3 C4 D" / "GWP- total kg CO2 eq. 2.363E+03 6.534E+01 ..."
    #   (the row label is split over lines in the text)                                          [VERIFIED]
    # Route, PDF PAGE 10: "TMT Bars at the Bhilai Steel Plant (BSP) are produced via the integrated Blast
    #   Furnace-Basic Oxygen Furnace (BF-BOF) route."
    # Period, PDF PAGE 7: "production of one metric tonne of TMT Bars at Bhilai Steel Plant for FY 2024-25".
    # Note: 2.363 kg/kg here and the gazette intensity of Tata Steel Jamshedpur (2.3630 t/t) are two different,
    #   unrelated numbers that happen to print the same digits. Both were read in their own source.
    ("SAIL_Bhilai_TMT_EPD-IES-0031470", "ore-based", "Steel Authority of India Limited (SAIL)", "TMT bars, Bhilai Steel Plant",
     "EPD-IES-0031470", Decimal("2363") / 1000, False, "SAIL Bhilai",
     "integrated Blast Furnace-Basic Oxygen Furnace (BF-BOF) route (EPD wording)", "FY 2024-25",
     "EPD_SAIL_Bhilai_TMT_IES-0031470.txt", "14"),
    # --- ARS Steels, EPD-IES-0012057 (PDF PAGE 1: "EPD registration number: EPD-IES-0012057") -----------
    # PDF PAGE 14, "Impact category Unit A1 A2 A3 A1-A3 C1 C2 C3 C4 D" /
    #   "GWP - total1) kg CO2e 1.01E+02 4.24E+01 4.48E+02 5.92E+02 0.00E+00 ..." -> printed A1-A3 total 592 kg/t
    #   (A1 + A2 + A3 from the rounded columns = 591.4)                                          [VERIFIED]
    # Route, PDF PAGE 11: "Electric induction furnace is used for melting the raw"; PDF PAGE 1: "(98% Recycled
    #   Steel)"; PDF PAGE 13: "Recycled Scrap -Local 816.6", "Recycled Scrap-Import 245.3", "Sponge Iron 14.2"
    #   (kg in a total of 1089).
    # Period, PDF PAGE 10: "Time representativeness: Jan 2024 to Dec 2024".
    # Scrap rule differs from the other EPDs, PDF PAGE 12: "100% of unsorted scrap is considered to carry an
    #   environmental burden".
    ("ARS_TMT_EPD-IES-0012057", "scrap", "ARS Steels and Alloy International Private Ltd", "TMT bars ARS 550D, ARS CRS 550D",
     "EPD-IES-0012057", Decimal("592") / 1000, False, None,
     "scrap in induction furnace; scrap carries a burden", "January - December 2024",
     "EPD_ARS_TMT_IES-0012057.txt", "14"),
    # --- Prakash Ferrous Industries, EPD-IES-0024699 (PDF PAGE 1: "EPD-IES-0024699") ---------------------
    # PDF PAGE 10, "Results for 1 Tonne of TMT Bars", "Impact indicators Unit A1-A3 C1 C2 C3 C4 D" /
    #   "GWP-total kg CO2 eq. 9.92E+02 0.00E+00 4.66E+00 ..."                                     [VERIFIED]
    # Route, PDF PAGE 4: "Steel scrap, high-grade sponge iron and other raw materials are melted in induction
    #   furnaces to produce billets." PDF PAGE 7: "The percentage of scrap here as an input is 107%." The EPD
    #   gives no mass of sponge iron, so the scrap share of the charge cannot be checked.
    # Period, PDF PAGE 7: "Time representativeness: April 2024- March 2025".
    ("Prakash_Ferrous_TMT_EPD-IES-0024699", "scrap", "Prakash Ferrous Industries Pvt. Ltd.", "TMT bars (500D, 550D, 550, 500D-CRS)",
     "EPD-IES-0024699", Decimal("992") / 1000, False, None,
     "scrap and sponge iron in induction furnaces; sponge iron amount not stated", "April 2024 - March 2025",
     "EPD_PrakashFerrous_TMT_IES-0024699.txt", "10"),
    # --- Jindal Steel & Power, EPD-IES-0008168 (PDF PAGE 1: "EPD registration number: EPD-IES-0008168") --
    # PDF PAGE 12, "Results per one metric ton of TMT Rebar", "Indicator Unit A1-A3 C1 C2 C3 C4 D" /
    #   "GWP- total kg CO2 eq. 4.45E+02 0.00E+00 7.62E+00 ..." (row label split over lines)        [VERIFIED]
    # Route, PDF PAGE 4: "100% recycled steel; 46% Pre-consumer scrap sourced internally & 54% post -consumer
    #   scrap ..., produced in an Electric Arc Furnace powered by a mix of electricity from 90% renewable energy
    #   and 10% from the grid." Angul plant (PDF PAGE 5).
    # Period, PDF PAGE 6: "The collection of the foreground data refers to the year 2023-24."
    # SPECIAL CASE: 90% renewable electricity is not typical of Indian scrap-based rebar. In the main case it
    #   is in the scrap group; one sensitivity row leaves it out.                                [ASSUMPTION]
    ("Jindal_TMT_100pct_scrap_EPD-IES-0008168", "scrap", "Jindal Steel & Power", "TMT rebar (100% recycled scrap)",
     "EPD-IES-0008168", Decimal("445") / 1000, True, None,
     "electric arc furnace, 100% scrap, 90% renewable electricity", "2023-24",
     "EPD_Jindal_TMT_100pct_scrap_IES-0008168.txt", "12"),
]
K_KEY, K_ROUTE, K_COMPANY, K_PRODUCT, K_REG, K_VAL, K_SPECIAL, K_PLANT, K_ROUTE_TXT, K_PERIOD, K_FILE, K_PAGE = range(12)
ROUTES = ["ore-based", "scrap"]
PLANT_ORDER = ["Tata Steel Jamshedpur", "JSW Steel Vijayanagar", "JSW Steel Dolvi", "SAIL Bhilai"]

# unit check: a value typed per tonne instead of per kg (or the reverse) stops the run
for e in EPD:
    assert Decimal("0.2") < e[K_VAL] < Decimal("5"), f"EPD value out of range (unit error?): {e[K_KEY]} = {e[K_VAL]}"
    assert e[K_ROUTE] in ROUTES, f"unknown route group: {e[K_KEY]}"
assert len({e[K_KEY] for e in EPD}) == len(EPD), "an EPD key is entered twice"
N_DOCS = len({e[K_REG] for e in EPD}); N_PRODUCERS = len({e[K_COMPANY] for e in EPD})
PLANT_EPD = {e[K_PLANT]: e for e in EPD if e[K_PLANT]}
for p in PLANT_ORDER:
    if p not in PLANT_EPD:
        raise KeyError(f"plant without an EPD value: {p}")
    assert PLANT_EPD[p][K_ROUTE] == "ore-based"

def obs_of(route, rel_sd=REL_SD_EPD, drop=None):
    return [(e[K_KEY], e[K_VAL], rel_sd * e[K_VAL]) for e in EPD if e[K_ROUTE] == route and e[K_KEY] != drop]

# =============================================================================
# PART 3 -- MAIN CASE: ONE UPDATE PER ROUTE GROUP (with the prior), AND WITHOUT THE PRIOR
# =============================================================================
MAIN = {r: pooled(MU0, SIGMA0, obs_of(r)) for r in ROUTES}
NOPRIOR = {r: pooled(None, None, obs_of(r)) for r in ROUTES}
SPREAD = {r: cochran(obs_of(r)) for r in ROUTES}              # spread test inside each group (EPDs only)
RANDOM = {r: random_effects(obs_of(r)) for r in ROUTES}       # random-effects estimate (EPDs only)   [ASSUMPTION]
REJECT = {r: SPREAD[r]["p"] < 0.05 for r in ROUTES}           # one common mean rejected at the 5% level

print("=" * 96)
print("STEP II -- route estimates for Indian rebar (kg CO2e / kg rebar, A1-A3, GWP-total)")
print("=" * 96)
print(f"Inputs: {len(EPD)} values, {N_DOCS} EPD documents, {N_PRODUCERS} producers. Prior {MU0} (SD {r3(SIGMA0)}); EPD SD {int(REL_SD_EPD*100)}% of each value.")
for r in ROUTES:
    vals = [e[K_VAL] for e in EPD if e[K_ROUTE] == r]
    for lab, R in (("with prior", MAIN[r]), ("no prior", NOPRIOR[r])):
        inside = min(vals) <= R["mu"] <= max(vals)
        print(f"  {r:10s} {lab:10s} n={len(vals)}  mu={r3(R['mu'])}  sd={r3(R['sd'])}  90%=[{r3(R['lo'])}, {r3(R['hi'])}]"
              f"  EPD range [{r3(min(vals))}, {r3(max(vals))}]  inside range: {inside}")
        # property checks: a weighted mean lies between the smallest and the largest input
        allx = R["x"]
        assert min(allx) <= R["mu"] <= max(allx)
        assert abs(sum(R["weights"]) - 1) < Decimal("1e-40")
    print(f"  {r:10s} prior weight {r1(MAIN[r]['weights'][0] * 100)}%")
print("\nSpread inside each group and random-effects estimate (DerSimonian-Laird, EPDs only, no prior):")
for r in ROUTES:
    H, E = SPREAD[r], RANDOM[r]
    print(f"  {r:10s} inputs [{r3(H['min'])}, {r3(H['max'])}]  Q={float(H['Q']):.3f} df={H['df']} p={H['p']:.3e} "
          f"I2={float(H['I2']):.3f}  common mean rejected at 5%: {'yes' if REJECT[r] else 'no'}  |  "
          f"random effects: tau2={float(E['tau2']):.6f} mu={r3(E['mu'])} sd={r3(E['sd'])} 90%=[{r3(E['lo'])}, {r3(E['hi'])}]")
for r in ROUTES:
    if REJECT[r]:
        print(f"  NOTE {r} group: one common mean is rejected. Its fixed-effect interval [{r3(MAIN[r]['lo'])}, "
              f"{r3(MAIN[r]['hi'])}] must NOT be quoted. Quote the input range [{r3(SPREAD[r]['min'])}, "
              f"{r3(SPREAD[r]['max'])}] and the random-effects interval [{r3(RANDOM[r]['lo'])}, {r3(RANDOM[r]['hi'])}].")

# Columns 1-6 are the fixed-effect estimate WITH the prior .
# The columns after them use the EPD inputs of the group only (no prior).
with open(f"{OUT}/route_estimates.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["route", "n", "mu_post", "sd_post", "ci90_lower", "ci90_upper",
                 "input_min", "input_max", "cochran_Q", "df", "p_value_chi2", "I2", "common_mean_rejected_at_5pct",
                 "tau2_dersimonian_laird", "mu_random_effects", "sd_random_effects",
                 "ci90_lower_random_effects", "ci90_upper_random_effects", "quote"])
    for r in ROUTES:
        R, H, E = MAIN[r], SPREAD[r], RANDOM[r]
        wr.writerow([r, len(obs_of(r)), f6(R["mu"]), f6(R["sd"]), f6(R["lo"]), f6(R["hi"]),
                     f6(H["min"]), f6(H["max"]), f6(H["Q"]), H["df"], f"{H['p']:.6e}", f6(H["I2"]),
                     "yes" if REJECT[r] else "no", f6(E["tau2"]), f6(E["mu"]), f6(E["sd"]), f6(E["lo"]), f6(E["hi"]),
                     ("input range and random-effects interval; NOT the fixed-effect interval" if REJECT[r]
                      else "fixed-effect interval may be quoted with the input range beside it")])

with open(f"{OUT}/epd_inputs.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["route_group", "input", "company", "product", "registration", "gwp_total_a1a3_kg_per_kg",
                 "sd_assumed", "weight_share_pct_in_route_main_case", "special_case", "plant", "data_period",
                 "source_file", "page", "route_note"])
    for r in ROUTES:
        R = MAIN[r]
        wr.writerow([r, "PRIOR (worldsteel global rebar, same for both routes)", "worldsteel", "global rebar", "",
                     f6(MU0), f6(SIGMA0), f6(R["weights"][0] * 100), "no", "", "2022 data release", "", "",
                     "global average over steelmaking routes; not route-specific"])
        for e in EPD:
            if e[K_ROUTE] != r:
                continue
            i = R["names"].index(e[K_KEY])
            wr.writerow([r, e[K_KEY], e[K_COMPANY], e[K_PRODUCT], e[K_REG], f6(e[K_VAL]), f6(R["s"][i]),
                         f6(R["weights"][i] * 100), "yes" if e[K_SPECIAL] else "no", e[K_PLANT] or "",
                         e[K_PERIOD], e[K_FILE], e[K_PAGE], e[K_ROUTE_TXT]])

# plant_epd_values.csv: for the paper's table that shows each plant's EPD value beside its GEI. The EPD value and
# the GEI are NEVER combined (no ratio, no scaling): they come from different accounting systems.
with open(f"{OUT}/plant_epd_values.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["plant", "epd_kgCO2e_per_kg", "epd_registration"])
    for p in PLANT_ORDER:
        wr.writerow([p, f6(PLANT_EPD[p][K_VAL]), PLANT_EPD[p][K_REG]])

# unrounded audit file
with open(f"{OUT}/audit_route_unrounded.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["route", "case", "source", "x", "sigma", "weight_1_over_sigma2", "w_times_x", "weight_share_pct"])
    for r in ROUTES:
        for case, R in (("main_with_prior", MAIN[r]), ("no_prior", NOPRIOR[r])):
            for n, x, s, w in zip(R["names"], R["x"], R["s"], R["w"]):
                wr.writerow([r, case, n, dec_str(x), dec_str(s), dec_str(w), dec_str(w * x), dec_str(w / R["W"] * 100)])
            for k in ("W", "S", "mu", "sd", "lo", "hi"):
                wr.writerow([r, case, k, dec_str(R[k]), "", "", "", ""])

# =============================================================================
# PART 4 -- SPREAD TEST (does one common value fit?)
# =============================================================================
ALL7 = [(e[K_KEY], e[K_VAL], REL_SD_EPD * e[K_VAL]) for e in EPD]
SPECIAL_KEYS = [e[K_KEY] for e in EPD if e[K_SPECIAL]]
HET = [("all seven values (both routes)", cochran(ALL7)),
       ("ore-based group", cochran(obs_of("ore-based"))),
       ("scrap group", cochran(obs_of("scrap"))),
       ("scrap group without the special case", cochran(obs_of("scrap", drop=SPECIAL_KEYS[0])))]
print("\nSpread test (Cochran Q; SD of each value = 10% of the value, an assumption):")
with open(f"{OUT}/heterogeneity.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["group", "k", "min_value", "max_value", "weighted_mean_no_prior", "cochran_Q", "df", "p_value_chi2",
                 "I2", "epd_rel_sd_assumed"])
    for lab, H in HET:
        wr.writerow([lab, H["k"], f6(H["min"]), f6(H["max"]), f6(H["mean_w"]), f6(H["Q"]), H["df"], f"{H['p']:.6e}",
                     f6(H["I2"]), f6(REL_SD_EPD)])
        print(f"  {lab:38s} k={H['k']}  Q={float(H['Q']):9.3f}  df={H['df']}  p={H['p']:.3e}  I2={float(H['I2']):.3f}"
              f"  range [{r3(H['min'])}, {r3(H['max'])}]")

# =============================================================================
# PART 5 -- SENSITIVITY PER ROUTE
# =============================================================================
sens = []     # route, case, n_epd, prior_used, R
for r in ROUTES:
    for rr in (Decimal("0.10"), Decimal("0.30"), Decimal("0.50")):                       # [ASSUMPTION]
        sens.append((r, f"prior_SD_{int(rr*100)}pct" + (" (main)" if rr == REL_SD0 else ""), "yes",
                     obs_of(r), pooled(MU0, rr * MU0, obs_of(r))))
    for rr in (Decimal("0.05"), Decimal("0.20")):                                        # [ASSUMPTION]
        sens.append((r, f"EPD_SD_{int(rr*100)}pct", "yes", obs_of(r, rr), pooled(MU0, SIGMA0, obs_of(r, rr))))
    # Equal absolute SD for all EPDs of the group (= REL_SD_EPD x mean EPD value of the group). The main case
    # uses an SD proportional to each value, which gives low values more weight; this case removes that. [ASSUMPTION]
    vals = [e[K_VAL] for e in EPD if e[K_ROUTE] == r]
    sd_eq = REL_SD_EPD * sum(vals) / len(vals)
    ob = [(n, v, sd_eq) for n, v, _ in obs_of(r)]
    sens.append((r, "EPDs_equal_SD", "yes", ob, pooled(MU0, SIGMA0, ob)))
    for e in EPD:
        if e[K_ROUTE] != r:
            continue
        ob = obs_of(r, drop=e[K_KEY])
        lab = f"without_{e[K_KEY]}" + (" (special case out)" if e[K_SPECIAL] else "")
        sens.append((r, lab, "yes", ob, pooled(MU0, SIGMA0, ob)))
    sens.append((r, "no_prior (EPDs only)", "no", obs_of(r), NOPRIOR[r]))
    ob = [(n, v, sd_eq) for n, v, _ in obs_of(r)]
    sens.append((r, "no_prior_EPDs_equal_SD (plain average)", "no", ob, pooled(None, None, ob)))
    if any(e[K_SPECIAL] for e in EPD if e[K_ROUTE] == r):
        ob = obs_of(r, drop=SPECIAL_KEYS[0])
        sens.append((r, "no_prior_without_special_case", "no", ob, pooled(None, None, ob)))

print("\nSensitivity per route:")
with open(f"{OUT}/route_sensitivity_cases.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["route", "case", "n_epd", "prior_used", "mu_post", "sd_post", "ci90_lower", "ci90_upper"])
    for r, lab, pu, ob, R in sens:
        wr.writerow([r, lab, len(ob), pu, f6(R["mu"]), f6(R["sd"]), f6(R["lo"]), f6(R["hi"])])
        print(f"  {r:10s} {lab:58s} mu={r3(R['mu'])} sd={r3(R['sd'])}")
# the leave-one-out rows that the paper names must be present (stops the run if an input key is renamed)
_labels = [lab for _, lab, _, _, _ in sens]
for _who in ("JSW_Dolvi", "Prakash_Ferrous", "Jindal_TMT"):
    _hit = [l for l in _labels if l.startswith("without_" + _who)]
    if len(_hit) != 1:
        raise KeyError(f"leave-one-out row 'without {_who}' not found exactly once in the sensitivity table")
    print(f"  leave-one-out row present: {_hit[0]}")
for r in ROUTES:
    mus = [R["mu"] for rr, _, _, _, R in sens if rr == r]
    print(f"  {r}: estimate ranges from {r3(min(mus))} to {r3(max(mus))} over {len(mus)} cases")

# =============================================================================
# PART 6 -- SCREENING TABLE: every Indian rebar-related EPD found, and why it was kept or left out.
# BEGIN SCREENING BLOCK. The four EPDs that were left out are named ONLY here and in epd_screening.csv.
# Their values are text for the table. They enter NO calculation in this or any later script.
# =============================================================================
SCREEN_DECISION = {"ore-based": ("ore-based", "rebar only; ore-based plant"),
                   "scrap": ("scrap", "rebar only; scrap melted in an induction furnace")}
screen_rows = []
for e in EPD:
    dec, why = SCREEN_DECISION[e[K_ROUTE]]
    if e[K_SPECIAL]:
        dec, why = "scrap special case", ("rebar only; 100% scrap in an electric arc furnace with 90% renewable electricity, "
                                          "not typical; in the scrap group in the main case, left out in one sensitivity row")
    if e[K_REG] == "S-P-04326":
        why += "; one EPD document gives two site values (Vijayanagar and Dolvi)"
    if e[K_KEY].startswith("JSW_Dolvi"):
        why = why.replace("ore-based plant", "ore-based plant (CONARC electric arc furnace, not BF-BOF)")
    if e[K_KEY].startswith("Prakash"):
        why += "; the EPD also names sponge iron as a raw material without an amount"
    screen_rows.append([e[K_COMPANY], e[K_PRODUCT], e[K_REG], e[K_ROUTE_TXT], e[K_PERIOD], f6(e[K_VAL]), dec, why,
                        e[K_FILE], e[K_PAGE]])
LEFT_OUT = [
    # Jindal Steel & Power, average long products. PDF PAGE 1: "S-P-09805 23-10-02 28-10-01 Global" and "Average
    #   Reinforcement Steel Bars, Round Bar, Wire Rod & RSC Steel Product". PDF PAGE 15, Table 12: "Climate Change -
    #   total [kg CO2 eq.] 3.45E+03 3.38E+01 ..." (columns A1-A3 A5 C1 C2 C3 C4 D). PDF PAGE 12: "the blast
    #   furnace/basic oxygen furnace route and the electric arc furnace route in Raigarh and Angul plants" and "The
    #   data was derived for the April 2021 - March 2022."
    ["Jindal Steel & Power", "Average of reinforcement steel bars, round bar, wire rod and RSC steel product", "S-P-09805",
     "blast furnace - basic oxygen furnace route and electric arc furnace route (billets from Raigarh and Angul)",
     "April 2021 - March 2022", "3.450000", "left out",
     "not rebar only: one average value for rebar, round bar, wire rod and another long product",
     "EPD_Jindal_avg_long_products_IES-0009805.txt", "15"],
    # Raajratna Metal Industries. PDF PAGE 1: "EPD REGISTRATION NUMBER S-P-09225". PDF PAGE 3: "Manufacturer and
    #   Exporter of Stainless Steel & Nickel Alloy Wires, Bright Bars, Reinforcement-bars" and "ASTM A 955/ SFS 1259".
    #   PDF PAGE 10: "Climate change - Total kg CO2 eq 2508.87 46.58 207.44 2762.88 ..." (A1, A2, A3, total).
    #   PDF PAGE 4: "refer to the financial year 2022-23".
    ["Raajratna Metal Industries Limited", "Steel rebar (stainless steel, ASTM A 955)", "S-P-09225",
     "hot rolling of purchased steel billet; steelmaking route not stated", "FY 2022-23", "2.762880", "left out",
     "stainless-steel rebar: a different product from carbon-steel rebar",
     "Raajratna_EPD_pdf.txt", "10"],
    # Peekay Rolling Mills. PDF PAGE 1: "EPD registration number: EPD- IES-0014568". PDF PAGE 11: "GWP- total kg CO2
    #   eq. 3.22E+03 0.00E+00 8.53E+00 ..." (A1-A3 first). PDF PAGE 10: "SPONGE IRON LUMPS & PELLETS 0.15083 MT",
    #   "MS BILLETS 0.3644 MT", "There is about  57.405% scrap input." PDF PAGE 5: scrap "is then melted by adding
    #   additives such as sponge iron" (furnace type not named). PDF PAGE 4: "Jan 2023- Dec 2023".
    ["Peekay Rolling Mills Pvt. Ltd.", "TMT bars", "EPD-IES-0014568",
     "scrap melted with sponge iron (furnace type not named), plus purchased steel billets", "January - December 2023",
     "3.220000", "left out",
     "route unclear: about 57% scrap input, sponge iron, and purchased billets of unstated origin in one value",
     "EPD_Peekay_TMT_IES-0014568.txt", "11"],
    # Bhuwalka Premier Group. PDF PAGE 1: "Registration number EPD-IES-0028749:003". PDF PAGE 20: "Climate change -
    #   total GWP-total kg CO2 eq. 1.70E+3 ND ND ..." (A1-A3 first). PDF PAGE 5: "The procured MS billets are cut to
    #   the required length and charged into a reheating furnace". PDF PAGE 8: "Scrap Steel 750", "Sponge Iron 250".
    #   PDF PAGE 10: "Raw material Database Ecoinvent v 3.12 2024 - 2025 Secondary data 0%" and total share of
    #   primary data "25.24%"; "The referance year is 2024 - 2025."
    ["Bhuwalka Premier Group", "TMT bars Fe 500, 500D, 550, 550D", "EPD-IES-0028749:003",
     "rolling mill only: purchased steel billets (declared content 75% scrap steel, 25% sponge iron) are reheated and rolled",
     "2024 - 2025", "1.700000", "left out",
     "steelmaking is not the producer's own: billets are bought and modelled with database data (primary data share 25.24%)",
     "EPD_Bhuwalka_TMT_IES-0028749.txt", "20"],
]
screen_rows += LEFT_OUT
with open(f"{OUT}/epd_screening.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["company", "product", "registration", "route_as_stated", "data_period", "gwp_total_a1a3_kg_per_kg",
                 "decision", "reason", "source_file", "page"])
    wr.writerows(screen_rows)
print(f"\nScreening table: {len(screen_rows)} rows ({len(EPD)} used, {len(LEFT_OUT)} left out).")
del LEFT_OUT, screen_rows
# END SCREENING BLOCK
# =============================================================================

# =============================================================================
# PART 7 -- FIGURE (1000 dpi + PDF TIFF, LZW)
# =============================================================================
def save(fig, path, dpi=1000):
    fig.tight_layout()
    fig.savefig(path, format="tiff", dpi=dpi, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(path.rsplit(".", 1)[0] + ".pdf", format="pdf", bbox_inches="tight")   # vector copy
    plt.close(fig)
    print("Wrote", path)

COL = {"ore-based": "#f16913", "scrap": "#41ab5d"}
fig, ax = plt.subplots(figsize=(8, 4.8), dpi=1000)
labels, y = [], len(EPD) - 1
for r in ROUTES:
    # the figure shows, for each group, the estimate that may be quoted (column "quote" of route_estimates.csv):
    # fixed effect with prior when one common mean is not rejected, random effects (EPDs only) when it is
    R = RANDOM[r] if REJECT[r] else MAIN[r]
    ax.axvspan(float(R["lo"]), float(R["hi"]), color=COL[r], alpha=0.18)
    ax.axvline(float(R["mu"]), color=COL[r], lw=1.6)
    for e in EPD:
        if e[K_ROUTE] != r:
            continue
        v, s = float(e[K_VAL]), float(REL_SD_EPD * e[K_VAL])
        ax.errorbar([v], [y], xerr=[s], fmt="D" if e[K_SPECIAL] else "o", color=COL[r], ecolor=COL[r], capsize=4, ms=6,
                    mec="black", mew=0.5)
        ax.text(v, y + 0.25, f"{v:.3f}", ha="center", fontsize=8)
        labels.append((y, e[K_KEY].replace("_", " ") + (" (special case)" if e[K_SPECIAL] else "")))
        y -= 1
ax.axvline(float(MU0), color="#666666", ls="--", lw=1.2)
ax.set_yticks([p for p, _ in labels]); ax.set_yticklabels([l for _, l in labels], fontsize=8)
ax.set_ylim(-0.7, len(EPD) - 0.3); ax.set_xlim(0, None)
ax.set_xlabel(r"kg CO$_2$e / kg rebar (A1-A3)", fontsize=12)
ax.set_title("Indian rebar EPD values and the two route estimates", fontsize=14, fontweight="bold", pad=15)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
def _fig_label(r, name):
    R = RANDOM[r] if REJECT[r] else MAIN[r]
    kind = "random effects, EPDs only" if REJECT[r] else "fixed effect, with prior"
    return f"{name} ({kind}): {float(R['mu']):.3f} [{float(R['lo']):.3f}, {float(R['hi']):.3f}]"
ax.legend(handles=[
    Line2D([], [], color=COL["ore-based"], marker="o", mec="black", mew=0.5, lw=1.6,
           label=_fig_label("ore-based", "Ore-based")),
    Line2D([], [], color=COL["scrap"], marker="o", mec="black", mew=0.5, lw=1.6,
           label=_fig_label("scrap", "Scrap")),
    Line2D([], [], color="#666666", ls="--", lw=1.2, label=f"Global prior, worldsteel ({float(MU0):.2f})"),
    Patch(fc="#999999", alpha=0.3, label=r"Shaded: 90% interval of the estimate named in the legend; bars: assumed $\pm 1\sigma$"),
], frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2)
save(fig, f"{OUT}/fig_route_estimates.tiff")

import platform, scipy
with open(f"{OUT}/software_versions_phase2_route_estimates_colab.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["package", "version"])
    wr.writerows([["python", platform.python_version()], ["scipy", scipy.__version__],
                  ["matplotlib", matplotlib.__version__]])

print("\nNo national value for Indian rebar is given: the coal-based sponge iron and induction furnace route has no usable EPD here.")
print("Step II done. Files in", OUT)
