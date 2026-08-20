function result=run_example()
%RUN_EXAMPLE Run the one self-contained baroclinic validation example.

setup_paths();
required={'sw_f','sw_dist','sw_pres','sw_bfrq','sw_pden'};
missing=required(cellfun(@(f) exist(f,'file')~=2,required));
if ~isempty(missing)
    error('Install the CSIRO Seawater Library first. Missing: %s', ...
        strjoin(missing,', '));
end
options=struct('save_output',false,'make_figure',false);
result=baroclinic_instability_58N22W(16,options);
fprintf('Example complete: maximum growth %.6g d^-1; 1/K scale %.6g km.\n', ...
    result.maximum_growth_per_day,result.inverse_wavenumber_scale_km);
end

