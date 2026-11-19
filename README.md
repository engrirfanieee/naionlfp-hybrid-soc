# Cross-chemistry SOC estimation for sodium-ion/LiFePO4 hybrid battery packs

Code and data for the article

> Habib MI. *Aging-aware cross-chemistry state-of-charge estimation for sodium-ion/LiFePO4 hybrid battery packs: a simulation analysis of the shared current-offset channel.* eTransportation (submitted), 2026.

The repository contains the physics-based (DFN) plant simulations, the identified equivalent-circuit models, all state estimators, the Monte Carlo and identifiability results, and the scripts that reproduce every table and figure of the article. Estimator names, case names and figure numbers are the same as in the article.

## Overview

A sodium-ion cell (Na3V2(PO4)2F3/hard carbon) and an LFP cell (LiFePO4/graphite) are connected in series and share one current. Both are simulated with Doyle–Fuller–Newman (DFN) models in PyBaMM, using the published parameter sets of Chayambuka et al. (2022) and Prada et al. (2012). The estimators use second-order equivalent-circuit models (ECMs) identified independently from pulse tests, so they never see the model that generated the data. Each estimator is evaluated over 20 Monte Carlo realisations of current offset, initial SOC errors and sensor noise in 11 test cases.

### Estimators (Section 3.9 of the article)

| Name | States | Shared offset *b* | Description |
|---|---|---|---|
| `CC` | – | – | Coulomb counting of the LFP cell, nominal capacity |
| `EKF-L` | zL, v1L, v2L | – | EKF on the LFP cell only |
| `IND-AA` | zN, v1N, v2N, zL, v1L, v2L, γN, γL, ρN, ρL | no | Two independent adaptive EKFs (capacity and resistance scales) |
| `IND-R` | zN, v1N, v2N, zL, v1L, v2L, ρN, ρL | no | Two independent resistance-adaptive EKFs, nominal capacities |
| `X-EKF` | zN, v1N, v2N, zL, v1L, v2L, b | yes | Coupled cross-chemistry EKF, nominal parameters |
| `X-EKF-R` | zN, v1N, v2N, zL, v1L, v2L, b, ρN, ρL | yes | Coupled, resistance-adaptive, nominal capacities |
| `AA-X-EKF` | all 11 states | yes | Proposed coupled, aging-aware EKF |

`IND-AA` and `IND-R` are implemented as one EKF over both cells, but they have no shared state. Their covariance therefore stays block-diagonal, and they are mathematically equivalent to two separate filters. All EKF variants use the same ECMs, measurement weighting and noise settings (`src/estimators.py`).

Ablations (Section 4.3): `X-EKF|b-known` and `IND-AA|b-known` (true offset subtracted, no offset state), `IND-AA-b` (each independent filter has its own offset state), `AA-X-EKF|b-calibrated` (offset constrained by a 1 mA calibration prior) and `AA-X-EKF|capacity-prior` (capacity scales initialised within ±3 % of the truth).

### Test cases (Table 2 of the article)

| Case | θ_LFP | θ_Na | ΔR_LFP | ΔR_Na | Profile | Plant files |
|---|---|---|---|---|---|---|
| A | 1.00 | 1.00 | 0 | 0 | UDDS | `plant_A_UDDS_*` |
| B | 0.95 | 0.88 | 0.15 | 0.30 | UDDS | `plant_B_UDDS_*` |
| C | 0.90 | 0.80 | 0.30 | 0.60 | UDDS | `plant_C_UDDS_*` |
| D | 0.85 | 0.70 | 0.45 | 0.90 | UDDS | `plant_D_UDDS_*` |
| E | 0.85 | 0.95 | 0.45 | 0.15 | UDDS | `plant_E_UDDS_*` |
| C-US06 | 0.90 | 0.80 | 0.30 | 0.60 | US06 | `plant_C_US06_*` |
| D-cap | 0.85 | 0.70 | 0 | 0 | UDDS | `plant_D_UDDS_*` |
| D-res | 1.00 | 1.00 | 0.45 | 0.90 | UDDS | `plant_A_UDDS_*` |
| F | D + inventory loss + halved kinetics | | 0.45 | 0.90 | UDDS | `plant_drv_F_*`, `fullcap_F_*` |
| A-inf | 1.00 | 1.00 | 0 | 0 | UDDS + voltage-limited C/4 charge + rest | `plant_inf2_A_*` |
| D-inf | 0.85 | 0.70 | 0.45 | 0.90 | UDDS + voltage-limited C/4 charge + rest | `plant_inf2_D_*` |

θ is the capacity retention (loss of active electrode area), and ΔR is the added series resistance as a fraction of the fresh ohmic resistance. ΔR is applied in post-processing: V = V_DFN − I·ΔR·R0,mean, with R0,mean = 39.5 mΩ (LFP) and 290 mΩ (Na-ion). In case F, the true capacities are the C/50 capacities of the degraded cells (1.819 Ah LFP, 1.471 Ah Na-ion). In the informative profiles, the charge is terminated when a terminal voltage (including ΔR) reaches 4.2 V (Na-ion) or 3.65 V (LFP). This took 34 min in A-inf and 17 min in D-inf; see `data/plant/inf2_A.json` and `inf2_D.json`.

## Repository structure

```
naion-lfp-hybrid-soc/
├── src/                         Python code (run every script from this folder)
│   ├── paths.py                 repository-relative paths
│   ├── plant.py                 DFN plant definition, aging, drive profiles
│   ├── run_plants.py            LFP plant runs; defines scenarios A–E (SCEN)
│   ├── run_na.py                sodium-ion plant runs (stepped solver)
│   ├── plant_ext.py             case F and the voltage-limited informative profiles (A-inf, D-inf)
│   ├── thermo_ocv.py            thermodynamic OCV from electrode potentials
│   ├── ident.py, ident2.py      pulse-test simulation and ECM identification
│   ├── modelsel.py              first- vs second-order ECM comparison (scenario A)
│   ├── estimators.py            all EKF variants and Coulomb counting
│   ├── evaluate.py              original six-estimator Monte Carlo comparison (Tables 6, 7, 9)
│   ├── evaluate_ext.py          all further Monte Carlo runs (Tables 5, 10–13)
│   ├── fim.py                   Fisher-information identifiability analysis (Table 4, Fig. 3)
│   ├── stats.py                 paired effect sizes, bootstrap CIs, Holm-adjusted tests (Table 8)
│   ├── exclusion_stats.py       regenerative-excursion statistics, exclusion-test effect sizes
│   ├── summarize.py             summary of the original comparison
│   ├── collect.py               all table values -> data/results/data_all.json
│   └── fig_*.py                 Figs. 1–6 and the graphical abstract
├── jobs/                        job lists for evaluate_ext.py
├── data/
│   ├── plant/                   DFN plant outputs (*.npz) and informative-profile metadata (*.json)
│   ├── plant_csv/               the same plant outputs as CSV
│   ├── model/                   identified OCV and ECM parameters (*.json)
│   └── results/
│       ├── montecarlo/          mc_<scenario>_<cycle>.json (original comparison)
│       ├── extended/            <case>__<method>__<tag>.json (all further runs)
│       ├── fim.json, effects.json, regen_exc.json, excl_eff.json
│       ├── summary.json         summary of the original comparison
│       └── data_all.json        every number used in the tables of the article
├── figures/                     Figs. 1–6 (PDF, PNG) and graphical abstract
├── requirements.txt, CITATION.cff, LICENSE
```

## Data description

### Plant data (`data/plant/*.npz`, `data/plant_csv/*.csv`)

Sampling interval 1 s. Columns (CSV) / keys (NPZ):

| CSV column | NPZ key | Unit | Meaning |
|---|---|---|---|
| `time_s` | `t` | s | time from start of test |
| `current_A` | `I` | A | applied current (positive = discharge) |
| `voltage_V` | `v` | V | DFN terminal voltage, **without** the added aging resistance ΔR |
| `discharge_capacity_Ah` | `q` | Ah | cumulative discharged charge |
| – | `npre` | samples | length of the preconditioning segment (7200 samples) |

Every test starts with a C/10 discharge that removes 15 % of nominal capacity, followed by a 30 min rest; the drive starts at sample 7200. The true SOC is z = 1 − q/Q, with Q = θ·Q_fresh (Q_fresh = 2.2958 Ah LFP, 2.3018 Ah Na-ion), except in case F (see above). `fullcap_F_<chem>` contains the C/50 discharge (q, V) and the capacity Q of the case-F cells.

### Model data (`data/model/`)

- `ecm2_params.json`: for each chemistry, the SOC grid, the thermodynamic OCV, the C/50 discharge and charge curves, the capacity `Q`, and the second-order pulse fits `ecm2`, one row per fit: [SOC, R0 (Ω), R1 (Ω), τ1 (s), R2 (Ω), τ2 (s), RMS residual (V), v1(0), v2(0)]. The estimators reject fits with an RMS residual above 25 mV.
- `ecm_params.json`: first-order fits (used only for the model-structure comparison).

### Results (`data/results/`)

- `montecarlo/mc_<scenario>_<cycle>.json`: for each estimator, a list of 20 realisations with `rmse`, `maxae`, `final` (LFP SOC, percentage points), `rmseN` (Na-ion SOC), `rdc_rmse`, `rdc_final` (mAh), the final capacity estimates `QN`, `QL` (Ah), the final offset estimate `b` (A) and the drawn `b_true`, `eN`, `eL`.
- `extended/<case>__<method>__<tag>.json`: the same per-realisation metrics plus `maxaeN`, `del_rmse` (RDC error against the charge actually delivered before the 2.0 V cutoff), `gap_mean`, `neesL`, `cov95L`, `neesN`, `cov95N` (consistency). In file names, `|` in the method name is written as `~`. Tags:
  - `base`: nominal settings
  - `excl`: measurement-exclusion test, Na-ion samples above 4.2 V discarded (Table 12)
  - `cons`: conservative covariance (σN = 75.6 mV, σL = 16 mV, SOC process noise 1e-8)
  - `sref0.5`, `sref2`, `kap10`, `kap1000`, `wide`, `sigN76`: tuning sensitivity (Table 13)
- `fim.json`: posterior standard deviations (square roots of the Cramér–Rao bound) and correlations for each profile, with the 11-parameter and the conditional 7-parameter results, and the 1 h window results.
- `effects.json`: paired differences against IND-AA and IND-R (Table 8).
- `regen_exc.json`, `excl_eff.json`: excursion statistics and exclusion-test effect sizes.
- `data_all.json`: every value used in the tables. Keys are `<case>|<method>|<tag>`, and each metric is [mean, std] over the 20 realisations.

## Reproducing the results

```bash
pip install -r requirements.txt
cd src

# 1) ECM identification (several minutes)            -> data/model/
python ident.py && python thermo_ocv.py && python ident2.py && python modelsel.py

# 2) DFN plant simulations (6–60 min each)            -> data/plant/
python run_plants.py A:UDDS B:UDDS C:UDDS D:UDDS E:UDDS C:US06
python run_na.py     A:UDDS B:UDDS C:UDDS D:UDDS E:UDDS C:US06
python plant_ext.py  F A D

# 3) Monte Carlo evaluation                           -> data/results/
for s in A B C D E; do python evaluate.py $s UDDS 20; done
python evaluate.py C US06 20
for j in ../jobs/*.json; do python evaluate_ext.py $j 2; done     # second argument: number of processes

# 4) identifiability, statistics, tables
python fim.py && python stats.py && python exclusion_stats.py && python summarize.py && python collect.py

# 5) figures                                          -> figures/
python fig_schem.py && python fig_ocv.py && python fig_identifiability.py && python fig_runs.py && python fig_montecarlo.py && python fig_ga.py
```

Steps 1–2 can be skipped, because the repository already contains their outputs. Steps 3–5 reproduce the stored results exactly. Each Monte Carlo job (20 realisations of one estimator in one case) takes about 1–4 min on one CPU core; the full set in step 3 takes several hours on a laptop. PyBaMM downloads the drive-cycle current profiles automatically (`pybamm.DataLoader`).

| Article item | Produced by | Data |
|---|---|---|
| Table 2 | `collect.py` (`_cases`) | `data_all.json` |
| Table 4, Fig. 3 | `fim.py`, `fig_identifiability.py` | `fim.json` |
| Table 5 (ablations) | `evaluate_ext.py jobs_ablations_and_X-EKF-R.json` | `extended/` |
| Tables 6, 7, 9, Fig. 6 | `evaluate.py`, `evaluate_ext.py`, `fig_montecarlo.py` | `montecarlo/`, `extended/` |
| Table 8 | `stats.py` | `effects.json` |
| Table 10 (capacities) | `evaluate_ext.py` | `extended/` (`QN`, `QL`) |
| Table 11 (aging modes, informative profile) | `jobs_aging_modes.json`, `jobs_informative_profile.json` | `extended/` |
| Table 12 (measurement exclusion) | `jobs_exclusion.json`, `exclusion_stats.py` | `extended/*__excl.json`, `excl_eff.json` |
| Table 13 (consistency, sensitivity) | `jobs_main_estimators.json`, `jobs_consistency_conservative.json`, `jobs_tuning_sensitivity.json` | `extended/` |
| Figs. 1, 2, 4, 5 | `fig_schem.py`, `fig_ocv.py`, `fig_runs.py` | `data/model/`, `data/plant/` |

## Random seeds

The same draws are used for every case and estimator:

- offsets *b* ~ U(±10 mA) and initial SOC errors ~ U(±15 percentage points): `numpy.random.default_rng(2026)`, 20 draws in the order (b, eN, eL);
- current and voltage noise of realisation *i* (*i* = 0…19): `default_rng(1000 + i)`;
- calibration-prior and capacity-prior draws of realisation *i*: `default_rng(1000 + i + 7919)`.

Sensor model: current noise σ = 10 mA, voltage noise σ = 1 mV, 1 mV quantisation. Accuracy metrics are evaluated from 1 h after the start of the drive.

## Software

Python 3.13, PyBaMM 26.10, NumPy 2.5, SciPy 1.18, Matplotlib 3.11 (see `requirements.txt`).

## Citation

Please cite the article above and this repository (see `CITATION.cff`).

## License

Code: MIT License. Data (`data/`) and figures (`figures/`): Creative Commons Attribution 4.0 International (CC BY 4.0). See `LICENSE`.

## References for the model parameters

- Chayambuka K, Jiang M, Mulder G, Danilov DL, Notten PHL. Physics-based modeling of sodium-ion batteries part I: experimental parameter determination. Electrochim Acta 2022;404:139726. https://doi.org/10.1016/j.electacta.2021.139726
- Chayambuka K, Mulder G, Danilov DL, Notten PHL. Physics-based modeling of sodium-ion batteries part II. Model and validation. Electrochim Acta 2022;404:139764. https://doi.org/10.1016/j.electacta.2021.139764
- Prada E, Di Domenico D, Creff Y, Bernard J, Sauvant-Moynot V, Huet F. Simplified electrochemical and thermal model of LiFePO4-graphite Li-ion batteries for fast charge applications. J Electrochem Soc 2012;159:A1508–19. https://doi.org/10.1149/2.064209jes
- Sulzer V, Marquis SG, Timms R, Robinson M, Chapman SJ. Python Battery Mathematical Modelling (PyBaMM). J Open Res Softw 2021;9:14. https://doi.org/10.5334/jors.309
