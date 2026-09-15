# Module 03: barotropic instability

This module evaluates the rotating shallow-water eigenvalue problem along the
complete 58 N section.

## Run

From the repository root:

```matlab
run('03_barotropic_instability/run_example.m')
```

## Inputs

- `data/reanaly_58_movavg.mat`: five daily velocity profiles from 9--13
  November 2020 on the complete 521-point, 0--520 km section.
- `data/smoothed_topography_on_xg.mat`: matching 521-point topography.

## Generated products

- `data/barotropic_instability_validation.mat`: modal diagnostics.
- `figure/barotropic_instability_validation.png`: three-panel diagnostic.
- `figure/barotropic_instability_validation.pdf`: vector figure.

Generated products are ignored by Git and can be recreated with
`run_example.m`.

## Configuration and scope

The example evaluates wavelengths from 10 to 300 km in 10-km increments,
uses a 1-km spatial grid, five-point topographic smoothing, Coriolis parameter
`1.23e-4 s^-1`, event threshold `0.1 d^-1`, and modal-period threshold of
20 days. Only the time dimension is reduced to five profiles, and the example
runs serially. 