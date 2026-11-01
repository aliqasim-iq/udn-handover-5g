# Supplementary material — "Decision Delay and Pre-Handover Service Gaps in Ultra-Dense 5G Networks: An Analytical Model and an Outage-Guarded A3 Trigger"

## tlmc/ — decision-delay model and trigger-level Monte-Carlo (TL-MC) simulator
- `model/udn_model.py` analytical model; `model/tl_mc.py` simulator; `model/run_all.py` experiments
  (`SEEDS=10 PROCS=4 python3 run_all.py kappa|val|scheme|sg|pp3|sens|lite`, needs numpy + scipy)
- `model/res_*.json` raw results used in Sections 4–6; `figs/*.py` regenerate Figures 1–5 (`tlmc/figs/fig1`–`fig5`)

## ns3_5g-lena/ — system-level evaluation (Section 7)
- Requires ns-3.48 + 5G-LENA v5.0 (commit 47a3adc). Copy the `.cc/.py` files into `scratch/`.
- `nr-udn-handover.cc`: scenario (square torus wrap-around, correlated LOS, OG-A3 guard, per-HO tracing)
- `nr-udn-handover-v2.cc`: identical, plus `--pktInterval` (experiment E5)
- `nr-udn-handover-v3.cc`: v2 plus the radio-level gap, MAC gap and outage-margin diagnostics, the A3 + return-block policy (`--policy=a3rb`) and the RLC-UM options (`--rlcReorderMs`, `--rlcOoo`); used for E7 and E8
- `nr-udn-rerun.py`: re-runs runs aborted by 5G-LENA (same seed truncated before the abort, or a new seed)
- `nr-udn-experiments.py`: run matrix E1–E4, E6 (OG-A3-lite), E7 (10 seeds x 60 s scheme comparison) and E8 (RLC t-Reordering 10 ms) (resumable, safe to run several instances)
  e.g. `python3 scratch/nr-udn-experiments.py --sets E4 --procs 4 --seeds 5 --ogThr -5.5 --ogCal -2.5`
- `run_e5.sh`: light-load experiment E5; `nr-udn-analysis.py`: statistics (`sys_stats.json`) and the system-level figures (`tlmc/figs/fig7`–`fig11`, see the mapping below)

## results/
- `5g-lena_raw/`: per-run summaries and per-HO records of all 5G-LENA runs
- `calibration_runs/`, `calib.json`: guard-threshold calibration (seeds 101–103, not used for evaluation)
- `sys_stats.json`: aggregated system-level statistics (mean, 95% CI)
- `run_logs/`: runner logs of E7/E8, including the aborted runs and how they were re-run

## Result tags
The TL-MC result files (`tlmc/model/res_*.json`) and the 5G-LENA run names use short tags:

| Tag | Article name |
|---|---|
| `fixed40`, `fixed100`, ... / `a3t40`, `a3t100`, ... | A3 with a fixed TTT of 40, 100, ... ms |
| `fixed100_nr` / `a3rbt40`, `a3rbt100` | A3 + RB |
| `mse` | MSE (oracle) |
| `ga` | MB-TTT, Eq. (19) |
| `og` / `og-5.5` | OG-A3 (guard threshold -5.5 dB) |
| `og_adapt` | OG-A3-adaptive, Eq. (20) |
| `og_nonr` | OG-A3-noRB |
| `og_lite` / `og-5.5-lite` | OG-A3-lite |
| `og_cho`, `og_ltm`, `ltm_fixed100` | OG-A3 + CHO-like, OG-A3 + LTM-like, LTM-like A3-100 |
| `og-2.5` | OG-A3 with the calibrated -2.5 dB threshold (E4 sensitivity point) |
| `-ro10` | RLC-UM t-Reordering 10 ms (E8) |

## Figures
`figures_paper/` holds the eight figures exactly as they appear in the article. The scripts write
their output to `tlmc/figs/` under the working names below.

| Article | Working file in `tlmc/figs/` | Script |
|---|---|---|
| Fig. 1 | `fig1_topology.png` | `fig_schematics.py` |
| Fig. 2 | `fig2_timeline.png` | `fig_schematics.py` |
| Fig. 3 | `fig3_validation.png` | `fig_model.py` |
| Fig. 4 | `fig4_kappa.png` | `fig_model.py` |
| Fig. 5 | `fig5_tradeoff.png` | `fig_schemes.py` |
| Fig. 6 | `fig10_sys_validation.png` | `nr-udn-analysis.py` |
| Fig. 7 | `fig11_sys_radio_kappa.png` | `nr-udn-analysis.py` |
| Fig. 8 | `fig7_sys_density.png` | `nr-udn-analysis.py` |

`fig6_speed.png`, `fig8_sys_speed_ttt.png` and `fig9_sys_schemes.png` are additional figures that
are not part of the article.

## License
- Python scripts, results and figures: MIT (see `LICENSE`).
- The ns-3 scenario files `ns3_5g-lena/*.cc` are derived from ns-3 / 5G-LENA and keep their
  `GPL-2.0-only` license, as stated in their headers.

## Citation
If you use this code or data, please cite the article
"Decision Delay and Pre-Handover Service Gaps in Ultra-Dense 5G Networks: An Analytical Model and an Outage-Guarded A3 Trigger" (submitted to *Computer Communications*).
