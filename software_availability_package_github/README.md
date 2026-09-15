# Iceland Current analysis software

This repository contains four MATLAB modules for the Iceland Current analysis.
Each module is self-contained: source code and compact validation inputs are
stored together, while generated MAT files and figures are excluded through
`.gitignore` and recreated by the module's `run_example.m`.

## Repository structure

```text
software_availability_package_github/
├── 01_data_assimilation/
│   ├── data/data_assimilation_validation.mat
│   ├── figure/
│   ├── fit_reanalysis_58N.m
│   ├── README.md
│   └── run_example.m
├── 02_eddy_detection_tracking/
│   ├── data/adt_validation_19930101_19930115.nc
│   ├── figure/
│   ├── eddy_detect.m
│   ├── first_birth.m
│   ├── README.md
│   └── run_example.m
├── 03_barotropic_instability/
│   ├── data/reanaly_58_movavg.mat
│   ├── data/smoothed_topography_on_xg.mat
│   ├── figure/
│   ├── time_var_reanalysis_58.m
│   ├── README.md
│   └── run_example.m
├── 04_baroclinic_instability/
│   ├── data/odv_58N_22p5W_validation.mat
│   ├── figure/
│   ├── baroclinic_instability_odv_58N22p5W.m
│   ├── README.md
│   └── run_example.m
├── scripts/
├── check_dependencies.m
├── run_validation.m
└── setup_paths.m
```

## Modules

| Module | Analysis | Bundled validation input | Generated products |
|---|---|---|---|
| 01 | Mooring-constrained 58 N velocity assimilation | 90-day GLORYS and four-mooring block | `reanaly_fit_58N.mat` and comparison PNG |
| 02 | ADT eddy detection and tolerance-LNN tracking | 15-day regional CMEMS/DUACS subset | detection/tracking MAT files and two PNGs |
| 03 | 58 N rotating shallow-water barotropic stability | five velocity profiles and matching topography | modal-diagnostic MAT, PNG, and PDF |
| 04 | ODV thermal-wind and continuously stratified QG stability | 38-level target profile and 3-by-3 sigma0 stencil | modal-diagnostic MAT, PNG, and PDF |

The compact inputs validate loading, numerical execution, saving, and plotting.
They do not replace the full observational and reanalysis archives used for
the manuscript.

## Setup and validation

Start MATLAB in the repository root:

```matlab
setup_paths
report = check_dependencies;
result = run_validation;
```

`run_validation` executes modules 01, 03, and 04 without writing generated
products. Module 02 is executed when OceanEddies and its required MathWorks
toolboxes are available; otherwise it is reported as skipped.

To generate the normal output files for one module:

```matlab
run('01_data_assimilation/run_example.m')
run('02_eddy_detection_tracking/run_example.m')
run('03_barotropic_instability/run_example.m')
run('04_baroclinic_instability/run_example.m')
```

## Dependencies

- MATLAB R2022b is the verified environment.
- Module 02 requires OceanEddies Eddyscan v2 and the Image Processing,
  Mapping, Statistics and Machine Learning, and Parallel Computing toolboxes.
- Module 03 can run serially; Parallel Computing Toolbox is optional.
- Module 04 uses `sw_f`, `sw_dist`, and `sw_pden` from the CSIRO MATLAB
  Seawater Library and uses `othercolor` for the supplied plotting style.

Place third-party MATLAB packages under `external/` or add an existing
installation to the MATLAB path. `setup_paths.m` automatically adds
`external/OceanEddies/` and `external/seawater/` when present.

## Validation-data builders

- `scripts/build_eddy_validation_data.py` rebuilds the module-02 ADT subset.
- `scripts/build_validation_data.py` rebuilds the two module-03 MAT inputs.
- `scripts/build_odv_baroclinic_validation_data.m` rebuilds the module-04 ODV
  MAT input from the original gridded binary fields.


## Repository hygiene

Only source code, documentation, compact validation inputs, and `.gitkeep`
placeholders belong in the repository. Generated MAT files, figures,
third-party packages, MATLAB autosave files, and operating-system metadata are
excluded by `.gitignore`.
