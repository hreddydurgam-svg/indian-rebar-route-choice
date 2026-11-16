# -*- coding: utf-8 -*-
"""
phase10_unit_checks_colab.py
============================
Two small checks on UNIT READINGS. It changes no earlier table and no earlier result.
It only reads tables written by phase 6 and phase 7 and writes one new table.

HOW TO RUN IN GOOGLE COLAB
  Run the five main scripts first (same runtime), or upload
     outputs_phase3/mass_energy_balance.csv
     outputs_phase4/twolimit_plant_options.csv
     outputs_phase4/twolimit_decisions.csv
  Then:  !python -X utf8 phase10_unit_checks_colab.py

CHECK 1 -- meaning of "ore share" in the capital cost
  Main model (phase 7): capex(s) = 184 + (1 - s) x 390 EUR per tonne, with s = scrap MASS share of the charge.
  Other reading: the electrolyser and the shaft furnace follow the ORE RATIO of the iron balance,
     ratio(s) = pellets(s) / pellets(0)   (phase 6, iron-balance scaling M),
     capex_alt(s) = 184 + ratio(s) x 390.
  The check solves all 100 limit cases again with capex_alt (same limits, same abatement, same scrap) and
  counts the cases in which the best plan changes.                                           [ASSUMPTION tested]

CHECK 2 -- "34.2 million tons" in the EY report
  Main model: read as 34.2 million metric tonnes.
  Other reading: short tons (1 short ton = 0.90718474 t). The check reports whether the main-case plan
  would still fit the scrap limit under that reading. The metric reading is the likely one for an Indian
  report; this check only shows how much rests on it.                                        [ASSUMPTION tested]

OUTPUT
  outputs_phase4/unit_checks.csv
"""
import csv
import itertools

OUT = "outputs_phase4"
CAPEX_EAF, CAPEX_ORE_PARTS = 184.0, 390.0        # EUR per tonne capacity; same values as phase 7 (Vogl et al., 2018)
SHORT_TON_T = 0.90718474                         # tonnes per short ton (exact by definition: 2000 lb x 0.45359237 kg)
TOL = 1e-9

opt = list(csv.DictReader(open(f"{OUT}/twolimit_plant_options.csv", encoding="utf-8")))
dec = list(csv.DictReader(open(f"{OUT}/twolimit_decisions.csv", encoding="utf-8")))
bal = list(csv.DictReader(open("outputs_phase3/mass_energy_balance.csv", encoding="utf-8")))

pel = {round(float(r["scrap_share"]), 2): float(r["value"]) for r in bal
       if r["stream"] == "iron ore pellets (iron-balance scaling M)"}
ratio = {s: pel[s] / pel[0.0] for s in (0.0, 0.50, 0.95)}
capex_main = {s: CAPEX_EAF + (1.0 - s) * CAPEX_ORE_PARTS for s in ratio}
capex_alt = {s: CAPEX_EAF + ratio[s] * CAPEX_ORE_PARTS for s in ratio}
CODE = {"N": "N", "S0": "0", "S50": "5", "S95": "9"}


def tables(elec, alt):
    t = {}
    for r in opt:
        if r["electricity"] != elec:
            continue
        o = r["option"]
        if o == "N":
            a = c = u = 0.0
        else:
            s = round(float(r["scrap_share"]), 2)
            a, u = float(r["abatement_Mt_per_yr"]), float(r["scrap_Mt_per_yr"])
            c = float(r["output_t_per_yr"]) * capex_alt[s] / 1e6 if alt else float(r["capital_million_EUR"])
        t.setdefault(r["plant_code"], {})[CODE[o]] = (a, c, u)
    return t


def best(elec, B, S, alt):
    t, top = tables(elec, alt), None
    for plan in itertools.product("N059", repeat=4):
        a = sum(t[p][o][0] for p, o in zip("JVDB", plan))
        c = sum(t[p][o][1] for p, o in zip("JVDB", plan))
        u = sum(t[p][o][2] for p, o in zip("JVDB", plan))
        if c <= B + TOL and u <= S + TOL:
            key = (-round(a, 9), round(c, 9), round(u, 9))      # tie rule of phase 7
            if top is None or key < top[0]:
                top = (key, "".join(plan), a)
    return top[1], top[2]


rows, n_changed = [], 0
for d in dec:
    B, S = float(d["capital_limit_million_EUR"]), float(d["scrap_limit_Mt"])
    p_main, a_main = best(d["electricity"], B, S, False)
    if p_main != d["best_plan"]:
        raise ValueError(f"this script does not reproduce the phase 7 plan in case {d['electricity']}, {d['f_B']}, {d['f_S']}")
    p_alt, a_alt = best(d["electricity"], B, S, True)
    if p_alt != p_main:
        n_changed += 1
        rows.append(["capex reading", f"{d['electricity']}, f_B={d['f_B']}, f_S={d['f_S']}", "best plan changes",
                     p_main, p_alt, f"{a_main:.6f}", f"{a_alt:.6f}", "Mt CO2e/yr"])

main = [d for d in dec if d["main_case"] == "1"][0]
p_main_alt, _ = best(main["electricity"], float(main["capital_limit_million_EUR"]), float(main["scrap_limit_Mt"]), True)
s_lim, s_used = float(main["scrap_limit_Mt"]), float(main["scrap_used_Mt"])
s_lim_short = s_lim * SHORT_TON_T

with open(f"{OUT}/unit_checks.csv", "w", newline="", encoding="utf-8") as f:
    wr = csv.writer(f)
    wr.writerow(["check", "item", "result", "main_reading", "other_reading", "value_main", "value_other", "unit"])
    for s in (0.0, 0.50, 0.95):
        wr.writerow(["capex reading", f"scrap share {s:.2f}", "ore share used for electrolyser and shaft",
                     "1 - s", "pellets(s) / pellets(0)", f"{1.0 - s:.6f}", f"{ratio[s]:.6f}", "-"])
        wr.writerow(["capex reading", f"scrap share {s:.2f}", "capital cost", "184 + (1 - s) x 390",
                     "184 + ratio x 390", f"{capex_main[s]:.6f}", f"{capex_alt[s]:.6f}", "EUR per tonne capacity"])
    wr.writerow(["capex reading", "all limit cases", "number of cases", "", "", str(len(dec)), "", "cases"])
    wr.writerow(["capex reading", "all limit cases", "cases in which the best plan changes", "", "",
                 "0", str(n_changed), "cases"])
    wr.writerow(["capex reading", "main case", "best plan", main["best_plan"], p_main_alt, "", "", "plan code"])
    wr.writerows(rows)
    wr.writerow(["scrap unit reading", "main case", "scrap limit", "million metric tonnes", "million short tons",
                 f"{s_lim:.6f}", f"{s_lim_short:.6f}", "Mt"])
    wr.writerow(["scrap unit reading", "main case", "scrap used by the best plan / scrap limit", "", "",
                 f"{s_used / s_lim:.6f}", f"{s_used / s_lim_short:.6f}", "-"])
    wr.writerow(["scrap unit reading", "main case", "best plan fits the scrap limit", "", "",
                 "yes" if s_used <= s_lim + TOL else "no", "yes" if s_used <= s_lim_short + TOL else "no", ""])

print("CHECK 1: ore ratio", {k: round(v, 6) for k, v in ratio.items()})
print("         capital cost, main", capex_main, "other reading", {k: round(v, 3) for k, v in capex_alt.items()})
print(f"         best plan changes in {n_changed} of {len(dec)} limit cases; main case: {main['best_plan']} -> {p_main_alt}")
print(f"CHECK 2: scrap used / limit = {s_used / s_lim:.4f} (metric), {s_used / s_lim_short:.4f} (short tons)")
print("Wrote", f"{OUT}/unit_checks.csv")
