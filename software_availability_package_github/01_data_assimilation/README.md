# Module 01: data assimilation

This module reconstructs the 58 N velocity section by combining a GLORYS
background with four aligned mooring observations.

## Run

From the repository root:

```matlab
run('01_data_assimilation/run_example.m')
```

All paths are resolved from the module directory.

## Input

`data/data_assimilation_validation.mat` contains:

- 90 daily samples from 1 November 2020 through 29 January 2021;
- 144 GLORYS section points;
- observations from UMM2, IB3, IB4, and IB5;
- nominal 0--215 m depth averages and centred 14-day observation means.

## Generated products

- `data/reanaly_fit_58N.mat`: assimilated 1-km velocity and 30-day background.
- `figure/reanalysis_fit_58N.png`: background, assimilated velocity, and
  observations for the selected validation date.

Generated products are ignored by Git and can be recreated with
`run_example.m`.

