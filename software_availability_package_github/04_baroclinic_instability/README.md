# Baroclinic-instability analysis at 58 N, 22 W

This module packages the author's original CMEMS/QG workflow as one MATLAB
function with paths resolved relative to the module directory.

## Run

Install the CSIRO MATLAB Seawater Library, then run from the repository root:

```matlab
setup_paths
result = baroclinic_instability_58N22W();
```

The default evaluates time index 16. To calculate without saving files or
opening a figure:

```matlab
options = struct('save_output',false,'make_figure',false);
result = baroclinic_instability_58N22W(16,options);
```

Default outputs are written to this module's `outputs/` directory.

The public function uses the local 58 N Coriolis and beta values, applies the
Bretherton/Smith boundary correction once, and averages bottom stratification
over the deepest finite positive N2 values rather than fill values below the
local seafloor.

