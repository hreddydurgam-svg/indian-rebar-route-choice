# Supplementary Information

**Route-specific carbon baselines and plant choice under capital and scrap limits for decarbonizing Indian steel rebar**

Himasai Kiran Reddy Durgam and Prasanna Venkatesan Ramani  
School of Civil Engineering, Vellore Institute of Technology, Vellore, India

Scripts and all output files: https://doi.org/10.5281/zenodo.24286235 (all versions); GitHub: https://github.com/hreddydurgam-svg/indian-rebar-route-choice

Scripts (run in this order): phase2_route_estimates_colab.py, phase3_route_model_colab.py, phase6_mass_energy_balance_colab.py, phase7_two_limit_selection_colab.py, phase8_twolimit_qaoa_colab.py, phase9_twolimit_figure_colab.py, phase10_unit_checks_colab.py; RUN_ORDER.txt.

Table S1. Screening of all Indian rebar EPDs found, with the reason each was kept or left out (epd_screening.csv). This is the only place where the four left-out EPDs are named.

Step II tables (outputs_phase2): epd_inputs.csv, route_estimates.csv, heterogeneity.csv, route_sensitivity_cases.csv, audit_route_unrounded.csv, plant_epd_values.csv. Figure: fig_route_estimates (Fig. 2 of the article).

Step III tables (outputs_phase3): plant_inputs.csv, scenario_table.csv, sweep_scrap.csv, sweep_grid.csv, plant_breakeven.csv, route_sensitivity.csv, route_vs_epd_evidence.csv, mass_energy_balance.csv, mass_energy_balance_sweep.csv, balance_checks.csv. Figures: fig_scrap_sweep (Fig. 3), fig_grid_sweep (Figure S1).

Step IV tables (outputs_phase4): twolimit_plant_options.csv, twolimit_decisions.csv, twolimit_summary.csv, twolimit_uncertainty_plans.csv, twolimit_uncertainty_summary.csv, twolimit_flip_thresholds.csv, twolimit_qaoa_main.csv, twolimit_qaoa_sweep.csv, twolimit_qaoa_summary.csv, unit_checks.csv. Figures: fig_twolimit_plans_a (Fig. 4), fig_twolimit_plans_b (Fig. 5). Five software version files (software_versions_*.csv). Figures 2 to 5 and Figure S1 were drawn again in Colab on 10 November 2026 at 1000 dpi, with vector PDF copies (Python 3.13.18, NumPy 2.4.6, SciPy 1.16.3, Matplotlib 3.10.0). This run wrote all phase 2 and phase 3 tables again; each is byte-identical to the run of 11 October 2026. Only the version files of phase 2 and phase 3 changed. The label of the 0.710 row in route_sensitivity.csv says UNVERIFIED; CEA version 22.0 (Appendix C, Table B, printed p. 28) later confirmed 0.710 for FY 2024-25 and 0.727 for FY 2023-24.

**Table S2. Mass and energy balance and route intensity per tonne of liquid steel (mass_energy_balance.csv, scenario_table.csv).**

| Item | 0% scrap | 50% scrap | 95% scrap |
|---|---|---|---|
| Pellets (kg) | 1,505 | 738.3 | 72.6 |
| Scrap (kg) | 0 | 536.2 | 1,001.9 |
| Hydrogen (kg) | 51.0 | 25.0 | 2.5 |
| Iron from pellets / scrap (kg) | 1,000 / 0 | 490.6 / 509.4 | 48.2 / 951.8 |
| Electrolyser (MWh) | 2.362 | 1.159 | 0.114 |
| EAF (MWh) | 0.753 | 0.710 | 0.671 |
| Heating remainder (MWh) | 0.365 | 0.179 | 0.018 |
| Electricity, sum of parts (MWh) | 3.480 | 2.048 | 0.803 |
| Electricity, straight line (MWh) | 3.480 | 2.074 | 0.808 |
| Route intensity, Indian grid (t CO2) | 2.808 | 1.786 | 0.865 |
| Route intensity, renewable (t CO2) | 0.176 | 0.176 | 0.176 |
| Capital cost (EUR/t capacity) | 574 | 379 | 203.5 |

Lime 50 kg and electrodes 2 kg at every scrap share. The heating remainder is scaled with the ore share (assumption). Pelletizing and rolling are not in the electricity rows.

QAOA settings (phase8_twolimit_qaoa_colab.py). Seed 12345. At each depth the angles start from a linear ramp (γ from 0.3 to 3.0, β from 0.75 to 0.075), from the best angles of the previous depth, and from 4 random starts. In the all-case run, the Grover-mixer method ranked a best-abatement plan first in 98 of 98 cases; its most likely plan was also the plan chosen by the exact solver in 68, and 31 cases have tied best plans, all under renewable electricity (twolimit_qaoa_summary.csv, twolimit_summary.csv).

**Table S3. Best plan in selected limit cases (plan code in the order Jamshedpur, Vijayanagar, Dolvi, Bhilai; N = not converted; 0, 5, 9 = converted at 0, 50, 95% scrap). Source: twolimit_decisions.csv.**

| Electricity | f_B | f_S | Best plan | Abatement (Mt/yr) | Capital (M EUR) | Scrap (Mt/yr) | Binding | Greedy plan (shortfall) |
|---|---|---|---|---|---|---|---|---|
| Indian grid | 0.50 | 0.10 | NNN5 | 7.74 | 2,151 | 3.04 | scrap | NNN5 (0%) |
| Indian grid (main) | 0.50 | 0.25 | NN55 | 16.08 | 5,526 | 7.82 | scrap | NNN9 (19.4%) |
| Indian grid | 0.50 | 0.50 | N559 | 31.07 | 8,979 | 16.76 | scrap | NN99 (5.1%) |
| Indian grid | 0.50 | 0.75 | 9N99 | 45.52 | 5,145 | 25.33 | scrap | 9N99 (0%) |
| Indian grid | 0.50 | 1.00 | 9959 | 57.90 | 9,097 | 32.95 | scrap | 5999 (2.9%) |
| Indian grid | 0.10 | 0.50 | NN9N | 16.53 | 1,812 | 8.92 | capital | NNN9 (21.6%) |
| Indian grid | 0.20 | 0.50 | NN99 | 29.49 | 2,967 | 14.61 | both | NN99 (0%) |
| Renewable | 0.50 | 0.25 | 05NN | 52.06 | 10,593 | 6.30 | both | NN09 (24.1%) |

Unit note: the grid factor is in g CO2 per kWh in scenario_table.csv, sweep_grid.csv and plant_breakeven.csv, and in t CO2 per MWh in route_sensitivity.csv and twolimit_plant_options.csv.
