# Iceland Current analysis software

This repository contains the MATLAB code used for the data assimilation, eddy
detection and tracking, barotropic-instability, and baroclinic-instability
analyses in the Iceland Current manuscript.

The repository is organized as four analysis modules. Author-specific absolute
paths have been removed. Generated products are written beneath `outputs/`, and
third-party software is expected beneath `external/`.

## Repository map

```text
iceland_current_analysis/
├── 01_data_assimilation/
│   └── fit_reanalysis_58N.m
├── 02_eddy_detection_tracking/
│   ├── eddy_detect.m
│   └── first_birth.m
├── 03_barotropic_instability/
│   └── time_var_reanalysis_58.m
├── 04_baroclinic_instability/
│   ├── baroclinic_instability_58N22W.m
│   └── data/
│       ├── uvts.nc
│       └── glo_topo.nc
├── data/
│   ├── reanaly_58_movavg.mat
│   └── smoothed_topography_on_xg.mat
├── outputs/
├── check_dependencies.m
├── run_example.m
└── setup_paths.m
```

## Analysis modules

1. `01_data_assimilation/fit_reanalysis_58N.m` reconstructs the 58 N velocity
   section from GLORYS and four mooring records. It forms nominal 0--215 m
   averages, applies a centered 14-day mean to the ADCP observations, and uses
   `gamma=400`, `lambda_background=0.5`, and `lambda_smoothness=50`. The
   reconstructed field is converted to the 30-day background state used by
   module 03.
2. `02_eddy_detection_tracking/eddy_detect.m` calls OceanEddies Eddyscan v2;
   `first_birth.m` then calls the tolerance-based LNN tracker and applies the
   study's area, lifetime, residence, and source-region criteria.
3. `03_barotropic_instability/time_var_reanalysis_58.m` performs the audited
   58 N rotating shallow-water eigenvalue analysis. Figure-4 events require a
   growth rate greater than 0.1 d^-1 and an Earth-fixed modal period shorter
   than 20 days.
4. `04_baroclinic_instability/baroclinic_instability_58N22W.m` packages the
   author's original continuously stratified QG workflow as one function. Its
   default reproduces time index 16 and the original `k-l` range of
   -2e-5--2e-5 m^-1.

## Quick start

Start MATLAB in the repository root and run:

```matlab
setup_paths
report = check_dependencies;
```

After installing the CSIRO Seawater Library, the included baroclinic example
can be run without author-specific files:

```matlab
result = run_example;
```

This evaluates module 04 at the original time index 16 without writing output
or opening a figure. Module-specific input and run instructions are given in
[`docs/INPUT_DATA.md`](docs/INPUT_DATA.md) and the comments at the beginning of
each MATLAB file.

## Third-party software

The code was verified with MATLAB R2022b (9.13). NetCDF input, linear least
squares, and dense eigenvalue calculations are provided by MATLAB; a separate
NetCDF library and the Optimization Toolbox are not required.

### CSIRO MATLAB Seawater Library

Module 04 requires release 3.3.1 of the
[CSIRO MATLAB Seawater Library](https://www.cmar.csiro.au/datacentre/ext_docs/seawater.html),
specifically `sw_f`, `sw_dist`, `sw_pres`, `sw_bfrq`, and `sw_pden`. Extract it
under `external/seawater/`, or add another installation to the MATLAB path.

The calculation uses the historical EOS-80 `sw_*` formulation. The TEOS-10
GSW toolbox is not a drop-in replacement; substituting it would require a
separate scientific reformulation and validation.

### OceanEddies

Module 02 requires the public
[OceanEddies](https://github.com/jfaghm/OceanEddies) MATLAB package. Install it
at `external/OceanEddies/`:

```bash
mkdir -p external
git clone https://github.com/jfaghm/OceanEddies.git external/OceanEddies
```

The author's working copy is labelled revision `7c4f9a4`. Before archival
release, retain that source snapshot or record the full Git commit identifier;
do not rely only on a moving branch. This study uses Eddyscan v2 `scan_single`
and `tolerance_track_lnn`, not the Python/Cython MHA tracker.

OceanEddies requires these MathWorks add-ons for the workflow used here:

- [Image Processing Toolbox](https://www.mathworks.com/products/image-processing.html)
  for connected-component, morphology, and region-property operations.
- [Mapping Toolbox](https://www.mathworks.com/products/mapping.html) for
  geographic distances and degree/kilometre conversions.
- [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html)
  for nearest-neighbour and radius searches.
- [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html)
  for the parallel OceanEddies v2 scan and faster module-03 calculations.
  Module 03 itself falls back to serial execution when this toolbox is absent.

The verified R2022b environment contained Image Processing Toolbox 11.6,
Mapping Toolbox 5.4, Statistics and Machine Learning Toolbox 12.4, and Parallel
Computing Toolbox 7.7. These are verified versions, not asserted minimums.

See [`docs/THIRD_PARTY_SOFTWARE.md`](docs/THIRD_PARTY_SOFTWARE.md) for the
dependency-to-module mapping and verification commands.

## Included and external data

The repository includes the two analysis-ready module-03 MAT files and the
58 N, 22 W module-04 NetCDF subset. The raw GLORYS/mooring inputs needed to
rerun module 01 and the daily CMEMS ADT input needed for module 02 are not
included. Their exact expected filenames and locations are listed in
[`docs/INPUT_DATA.md`](docs/INPUT_DATA.md).

The included NetCDF files retain their CMEMS/GLORYS metadata. Authors should
confirm dataset redistribution and acknowledgement requirements before making
the GitHub repository public.

## GitHub upload

`04_baroclinic_instability/data/uvts.nc` is approximately 79 MB. It is below
GitHub's 100 MiB hard limit for regular Git objects, but above the browser
upload limit. Upload this repository with Git on the command line, not by
dragging the files into the GitHub web interface. Detailed commands are in
[`docs/UPLOAD_TO_GITHUB.md`](docs/UPLOAD_TO_GITHUB.md).

## Reproducibility and attribution

- Do not commit generated outputs, MATLAB autosave files, or operating-system
  metadata; these are excluded by `.gitignore`.
- Preserve the source identifiers and acknowledgement text for reused CMEMS,
  GLORYS, OSNAP, OceanEddies, and Seawater resources.
- Third-party code is linked rather than silently vendored. Preserve its
  original licence and attribution information.
- Before public release, the authors must select a licence for their own code
  and add the corresponding `LICENSE` file. No licence has been assumed here.

