# Module 04: baroclinic instability

This module packages the ODV thermal-wind and QG stability calculation from
the latter part of `2nd_round/icelandcurrent.m`. The compact input is labelled
22.5 W, 58 N and retains the source extraction indices `(76,279)`.

## Run

From the repository root:

```matlab
run('04_baroclinic_instability/run_example.m')
```

## Input

`data/odv_58N_22p5W_validation.mat` contains 38 temperature, salinity,
pressure, and depth levels plus the 3-by-3 sigma0 stencil used for the central
horizontal density gradients.

## Generated products

- `data/baroclinic_instability_odv_validation.mat`: density, thermal-wind,
  eigenvalue, eigenvector, amplitude, phase, and local-maximum results.
- `figure/baroclinic_instability_odv_validation.png`: growth-rate contour.
- `figure/baroclinic_instability_odv_validation.pdf`: vector figure.

Generated products are ignored by Git and can be recreated with
`run_example.m`.

## Source configuration

The module retains the supplied program's calculation sequence:

- potential density referenced to 40 dbar;
- zero relative velocity at level 22;
- thermal-wind gravity `9.8 m s^-2`;
- the first 31 levels in the QG calculation;
- reference density `1024 kg m^-3` and QG gravity `9.81 m s^-2`;
- 100 values each of `k` and `l` from `-2e-3` to `2e-3 m^-1`;
- the original local-maximum search and triangular plot markers.

The source pressure selection `p(1,279,:)` is retained. The data extraction
location and plotting layout are unchanged.

## Dependencies and data rebuild

The calculation requires `sw_pden`, `sw_f`, and `sw_dist` from the CSIRO
MATLAB Seawater Library and `othercolor` for plotting. Rebuild the compact
input entirely in MATLAB with:

```matlab
addpath('scripts')
build_odv_baroclinic_validation_data('/path/to/gridded')
```
