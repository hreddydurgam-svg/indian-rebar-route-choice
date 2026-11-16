# Supplementary code and data

Supplementary Information for the article
"Route-specific carbon baselines and plant choice under capital and scrap limits for decarbonizing Indian steel rebar"
(Himasai Kiran Reddy Durgam and Prasanna Venkatesan Ramani, School of Civil Engineering, Vellore Institute of Technology, India).

The folder holds the seven Python scripts of the study and every table and figure they wrote.
`SUPPLEMENTARY_INFORMATION.md` gives the list of supplementary items, Table S2, Table S3 and the QAOA settings.
All numbers in the article come from these output files or from the cited source documents.

## How to run

Run the scripts in Google Colab, in one runtime, in the order given in `RUN_ORDER.txt`.
No input file is needed. Each script reads only files written by the scripts before it.
Example: `!python -X utf8 phase2_route_estimates_colab.py`

| Step | Script | Writes to |
|---|---|---|
| II | `phase2_route_estimates_colab.py` | `outputs_phase2/` |
| III | `phase3_route_model_colab.py` | `outputs_phase3/` |
| III | `phase6_mass_energy_balance_colab.py` | `outputs_phase3/` |
| IV | `phase7_two_limit_selection_colab.py` | `outputs_phase4/` |
| IV | `phase8_twolimit_qaoa_colab.py` | `outputs_phase4/` |
| IV | `phase9_twolimit_figure_colab.py` (draws Figs. 4 and 5 only) | `outputs_phase4/` |
| IV | `phase10_unit_checks_colab.py` (tests two unit readings) | `outputs_phase4/` |

## Output files

**Step II, `outputs_phase2/`**
- `epd_screening.csv`: every Indian rebar EPD found, with the reason each was kept or left out (Table S1).
- `epd_inputs.csv`, `route_estimates.csv`, `heterogeneity.csv`, `route_sensitivity_cases.csv`, `audit_route_unrounded.csv`, `plant_epd_values.csv`.
- `fig_route_estimates.tiff` / `.pdf`: Fig. 2.

**Step III, `outputs_phase3/`**
- `plant_inputs.csv`, `scenario_table.csv`, `sweep_scrap.csv`, `sweep_grid.csv`, `plant_breakeven.csv`, `route_sensitivity.csv`, `route_vs_epd_evidence.csv`.
- `mass_energy_balance.csv`, `mass_energy_balance_sweep.csv`, `balance_checks.csv` (Table S2 uses `mass_energy_balance.csv` and `scenario_table.csv`).
- `fig_scrap_sweep.tiff` / `.pdf`: Fig. 3. `fig_grid_sweep.tiff` / `.pdf`: Figure S1.

**Step IV, `outputs_phase4/`**
- `twolimit_plant_options.csv`, `twolimit_decisions.csv`, `twolimit_summary.csv`, `twolimit_uncertainty_plans.csv`, `twolimit_uncertainty_summary.csv`, `twolimit_flip_thresholds.csv`.
- `twolimit_qaoa_main.csv`, `twolimit_qaoa_sweep.csv`, `twolimit_qaoa_summary.csv`.
- `unit_checks.csv`.
- `fig_twolimit_plans_a.tiff` / `.pdf`: Fig. 4. `fig_twolimit_plans_b.tiff` / `.pdf`: Fig. 5.

**Software versions:** five files `software_versions_*.csv`, one for each of phases 2, 3, 6, 7 and 8.

## Runs and software

- Phases 6, 7 and 8: Google Colab, 11 October 2026. Python 3.13.16 (phases 6, 7, 8); NumPy 2.1.3 (phases 7, 8); SciPy 1.16.3 (phase 8).
- Phases 2 and 3: run again in Google Colab on 10 November 2026 to draw the figures at 1000 dpi with vector PDF copies. Python 3.13.18 (phases 2, 3); SciPy 1.16.3 and Matplotlib 3.10.0 (phase 2); NumPy 2.4.6 and Matplotlib 3.10.0 (phase 3). All phase 2 and phase 3 tables of this run are byte-identical to those of the 11 October 2026 run.
- The versions above are those recorded in the five `software_versions_*.csv` files.
- Phase 9 (figures only) was run in Google Colab on 10 November 2026; phase 10 was run in Google Colab on 12 October 2026. These two scripts write no version file.
- The emulated QAOA (phase 8) is a classical state-vector emulation with NumPy and SciPy. No quantum hardware and no quantum software library was used.

## Notes

- **Comments edited after the runs.** Some code comments and docstrings were edited after the Colab runs: process notes were removed, and the note on the 0.710 grid factor was brought up to date. No code line was changed. For each script, the Python syntax tree with docstrings set aside is identical to that of the script that was run.
- **Grid factor 0.710.** The label of this row in `route_sensitivity.csv` says "UNVERIFIED", because the value had not been checked in a CEA document when the script was run. CEA User Guide Version 22.0 (Appendix C, Table B) later confirmed 0.710 t CO2/MWh for FY 2024-25 and 0.727 t CO2/MWh for FY 2023-24. The main value in all scripts is 0.727 (FY 2023-24).
- **worldsteel prior.** The prior of 1.73 kg CO2e/kg rebar is one aggregated result from the worldsteel LCI data (2022 release), used under worldsteel terms. The worldsteel data report and inventory are not shared here.
- **Sources.** EPD, CEA, CCTS (MoEFCC), EY and Ministry of Steel values come from the public documents cited in the article and in the script comments. The `source_file` columns in some tables name the text files that the authors extracted from those documents. These text files are not part of this package.

## License

- Scripts (`*.py`): MIT License, see `LICENSE`.
- Output tables and figures (`outputs_phase2/`, `outputs_phase3/`, `outputs_phase4/`): CC BY 4.0, see `LICENSE-DATA.md`.
