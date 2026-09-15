# Module 02: eddy detection and tracking

This module runs the manuscript's ADT eddy workflow: Eddyscan v2 detection,
tolerance-based LNN tracking, track filtering, and birth-location analysis.

## Run

From the repository root:

```matlab
run('02_eddy_detection_tracking/run_example.m')
```

## Input

`data/adt_validation_19930101_19930115.nc` contains 15 consecutive days of
CMEMS/DUACS absolute dynamic topography on the complete regional grid. Only
the time dimension is reduced.

## Generated products

- `data/eddy_detections_validation.mat`: daily Eddyscan v2 detections.
- `data/eddy_tracking_validation.mat`: accepted tracks and birth locations.
- `figure/eddy_detection_validation.png`: final-day detections.
- `figure/eddy_tracking_validation.png`: tracks and normalized birth counts.

Generated products are ignored by Git and can be recreated with
`run_example.m`.

## Validation settings

| Parameter | Compact example | Full workflow |
|---|---:|---:|
| Physical-area prefilter | none | 8,000 km^2 |
| Tracker gap tolerance | 1 day | 3 days |
| Minimum track lifetime | 4 days | 27 days |
| Minimum source-region residence | 1 day | 8 days |
| Birth histogram spacing | 0.5 degrees | 0.2 degrees |

Eddyscan retains the v2 algorithm and a nine-pixel minimum detection area.


## Dependencies

Module 02 requires OceanEddies plus MATLAB's Image Processing, Mapping,
Statistics and Machine Learning, and Parallel Computing toolboxes. Place
OceanEddies under `external/OceanEddies/` or add it to the MATLAB path.
