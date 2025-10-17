# slag-kinetic

Code and data for kinetic analysis of zinc oxide reduction in slag, supporting:

The catalytic effect of the metal bath on the zinc oxide reduction. Results in Engineering, 16 (2022) 100514. https://doi.org/10.1016/j.rineng.2022.100514

This repository contains a small Python library for:
- Parsing elemental analyses (e.g., SEM‑EDX) and computing phase compositions
- Assembling per‑trial datasets with time and mass information
- Converting analyses to mol courses via inert tracer normalization
- Fitting ODE‑based kinetic models (SciPy) for ZnO/FeO systems
- Producing example plots, residuals, and reaction‑rate vs. composition views

## Repository Structure
- `nemlib/`
  - `analysis.py` — Elemental/Phase analyses, file parsers, info merging
  - `chemical.py` — Element/molecule utilities (masses, parsing, conversions)
  - `model.py` — Kinetic models (first‑order, ZnO crucible, ZnO–FeO coupled)
  - `modelhelper.py` — Trial dimensions, tracer‑based mol‑course, assemblers
  - `trial.py` — Trial organization, selectors, time axis utilities
- `slag_trials_2021/`
  - `data/sem.xlsx` — Elemental analyses (wt.% input expected)
  - `data/time_table.xlsx` — Sampling times by sample ID
  - `data/trial_info.xlsx` — Per‑trial metadata (e.g., inert mass, setup)
  - `setup.py` — Helpers to assemble trials and produce figures
  - `plots/` — Example output (e.g., `detail_plots.svg`)
- `tests/` — Unit tests for parsing, conversions, and model behavior

## Requirements
- Python 3.9+
- Packages: `numpy`, `pandas`, `scipy`, `matplotlib`, `openpyxl`, `pytest`

Install:

```bash
pip install numpy pandas scipy matplotlib openpyxl pytest
```

## Quickstart
Minimal example using the included dataset and a ZnO crucible model:

```python
from nemlib.analysis import PhaseConstructor, ElementalAnalysesFile, InfoOrganizer
from nemlib.trial import SlagReductionTrial, VxPySelector
from nemlib.modelhelper import TrialAssembler, ModelCompounds, TrialDimensions
from nemlib.model import ZnOCCrucibleModel
import numpy as np

# 1) Load elemental -> derive phases
phases = ["CaF2","CaO","SiO2","Al2O3","MgO","Na2O","K2O","FeO","ZnO"]
pc = PhaseConstructor(phases=phases, ignore=["O","C"])
data = ElementalAnalysesFile("slag_trials_2021/data/sem.xlsx").phased(pc)

# Add per-sample time and per-trial info
data.add_info(InfoOrganizer("slag_trials_2021/data/time_table.xlsx", "sample", ["time"]))
trial_info = InfoOrganizer("slag_trials_2021/data/trial_info.xlsx", "trial",
                           ["inert mass","zno mass","metal mass","setup"])

# 2) Assemble a kinetic trial
assembler = TrialAssembler(
    trial_info=trial_info,
    inert_compounds=["CaO","SiO2","Al2O3","MgO","CaF2","K2O","Na2O"],
    dimensions=TrialDimensions,   # computes areas from crucible dia, mass, density
    slag_density=2.5,             # g/cm^3 (used for slag height)
    model=ZnOCCrucibleModel,
    model_compounds=ModelCompounds(["ZnO","FeO"]),
)

trial = SlagReductionTrial(VxPySelector("V29")(data))
trial = assembler.do(trial, exclude=[])     # drop sample indices here if needed
trial.model.shift_to_y0(12)                 # align to x_ZnO(t=0)=12 mol-%
t, zno_pct = trial.model_y(pct=True)        # fitted ZnO mol-% vs time

print("Fit report:", trial.model.fit_report)  # mean ± stdev of residuals
```

For figure generation and model comparisons, see `slag_trials_2021/setup.py`:
- `MetalBathSeries()` sets up common loaders and assemblers
- `SummaryPlot` and `DetailPlot` create overview and per‑trial plots (including deviations and reaction‑rate vs composition)

## Tests
Run unit tests:

```bash
pytest -q
```

## Data Notes
- Elemental inputs are interpreted as wt.% and converted internally.
- Inert tracer normalization uses the specified inert phases and the initial inert mass per trial.
- Example Excel inputs:
  - `sem.xlsx`: elemental analyses (samples in rows)
  - `time_table.xlsx`: columns `sample`, `time`
  - `trial_info.xlsx`: columns `trial`, `inert mass`, `zno mass`, `metal mass`, `setup`

## Citation
Please cite:

The catalytic effect of the metal bath on the zinc oxide reduction. Results in Engineering, 16 (2022) 100514. https://doi.org/10.1016/j.rineng.2022.100514

## License
- Code: MIT License (see `LICENSE`)
- Data: CC BY 4.0 recommended for datasets in `slag_trials_2021/data` (see `LICENSE-DATA`)
