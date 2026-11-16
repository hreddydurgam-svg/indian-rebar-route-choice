# -*- coding: utf-8 -*-
"""
phase9_twolimit_figure_colab.py
===============================
Draws TWO figures from tables that the earlier scripts wrote. It computes no new result.

HOW TO RUN IN GOOGLE COLAB
  Run phase7_two_limit_selection_colab.py and phase8_twolimit_qaoa_colab.py first (same runtime), or upload
     outputs_phase4/twolimit_decisions.csv
     outputs_phase4/twolimit_qaoa_main.csv
  Then:  !python -X utf8 phase9_twolimit_figure_colab.py

OUTPUT
  outputs_phase4/fig_twolimit_plans_a.tiff/.pdf  (1000 dpi, LZW; vector PDF)
     (a):  best plan and its abatement for every capital limit and scrap limit, India grid
            (plan code: plant order Jamshedpur, Vijayanagar, Dolvi, Bhilai; N = not converted,
             0 / 5 / 9 = converted at 0 / 50 / 95 % scrap); the main case is framed
     outputs_phase4/fig_twolimit_plans_b.tiff/.pdf
     (b): probability of the best plan against depth, main case, for the four emulated methods
"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "outputs_phase4"
DEC = list(csv.DictReader(open(f"{OUT}/twolimit_decisions.csv", encoding="utf-8")))
QA = list(csv.DictReader(open(f"{OUT}/twolimit_qaoa_main.csv", encoding="utf-8")))

rows = [r for r in DEC if r["electricity"] == "India grid"]
FB = sorted({r["f_B"] for r in rows}, key=float)
FS = sorted({r["f_S"] for r in rows}, key=float)
A = np.zeros((len(FB), len(FS)))
PLAN = {}
MAIN = None
for r in rows:
    i, j = FB.index(r["f_B"]), FS.index(r["f_S"])
    A[i, j] = float(r["abatement_Mt_per_yr"])
    PLAN[(i, j)] = r["best_plan"]
    if r["main_case"] == "1":
        MAIN = (i, j)

fig, ax = plt.subplots(figsize=(6.0, 4.8), dpi=1000)
figb, bx = plt.subplots(figsize=(5.2, 4.8), dpi=1000)
im = ax.imshow(A, origin="lower", aspect="auto", cmap="YlGn")
for (i, j), p in PLAN.items():
    ax.text(j, i, f"{p}\n{A[i, j]:.1f}", ha="center", va="center", fontsize=7,
            color="white" if A[i, j] > 0.6 * A.max() else "black")
if MAIN is not None:
    ax.add_patch(plt.Rectangle((MAIN[1] - 0.5, MAIN[0] - 0.5), 1, 1, fill=False, ec="#d7301f", lw=2.2))
ax.set_xticks(range(len(FS))); ax.set_xticklabels(FS)
ax.set_yticks(range(len(FB))); ax.set_yticklabels(FB)
ax.set_xlabel("Scrap limit (share of 34.2 Mt)", fontsize=11)
ax.set_ylabel("Capital limit (share of cost of all four plants at 0% scrap)", fontsize=9)
ax.set_title("Best plan and abatement, India grid", fontsize=11, fontweight="bold")
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
cb.set_label(r"Abatement (Mt CO$_2$e / yr)", fontsize=9)

STYLE = {"G": ("#238b45", "o", "Grover mixer on feasible plans (G)"),
         "O_hinge": ("#2171b5", "s", "Linear penalty, smallest separating (O_hinge)"),
         "O_fixed": ("#f16913", "^", "Squared penalty, fixed rule (O_fixed)"),
         "O_min": ("#6a51a3", "v", "Squared penalty, smallest separating (O_min)")}
for m, (c, mk, lab) in STYLE.items():
    pts = sorted((int(r["p"]), float(r["P_best"])) for r in QA if r["method"] == m)
    bx.plot([p for p, _ in pts], [v for _, v in pts], color=c, marker=mk, lw=1.6, ms=5, label=lab)
eq = [float(r["P_best"]) for r in QA if r["method"] == "equal_chance_feasible"][0]
bx.axhline(eq, color="#666666", ls="--", lw=1.1, label=f"Equal chance over feasible plans ({eq:.3f})")
bx.set_xlabel("QAOA depth p", fontsize=11); bx.set_ylabel("Probability of the best plan", fontsize=11)
bx.set_ylim(0, 1); bx.set_title("Main case, emulated QAOA", fontsize=11, fontweight="bold")
bx.spines["top"].set_visible(False); bx.spines["right"].set_visible(False)
bx.legend(frameon=False, fontsize=7.5, loc="upper left")
for f, tag in ((fig, "a"), (figb, "b")):
    f.tight_layout()
    f.savefig(f"{OUT}/fig_twolimit_plans_{tag}.tiff", format="tiff", dpi=1000, bbox_inches="tight",
              pil_kwargs={"compression": "tiff_lzw"})
    f.savefig(f"{OUT}/fig_twolimit_plans_{tag}.pdf", format="pdf", bbox_inches="tight")
    print("Wrote", f"{OUT}/fig_twolimit_plans_{tag}.tiff and .pdf")
