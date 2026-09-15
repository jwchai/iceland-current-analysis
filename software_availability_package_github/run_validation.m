function result=run_validation()
%RUN_VALIDATION Exercise all available modules with compact sample inputs.
%

root=setup_paths();
result=struct();

options=struct('make_figure',false,'save_output',false);
result.data_assimilation=fit_reanalysis_58N([],[],options);
fprintf('Module 01 validation complete: %d days assimilated.\n', ...
    numel(result.data_assimilation.time));

eddy_required={'scan_single','tolerance_track_lnn','regionprops', ...
    'distance','km2deg','rangesearch','knnsearch'};
eddy_missing=eddy_required(cellfun(@(f) exist(f,'file')~=2,eddy_required));
if isempty(eddy_missing)
    temporary_dir=tempname;
    mkdir(temporary_dir);
    cleanup_eddy=onCleanup(@() remove_validation_temp(temporary_dir)); %#ok<NASGU>
    input_file=fullfile(root,'02_eddy_detection_tracking','data', ...
        'adt_validation_19930101_19930115.nc');
    detection_file=fullfile(temporary_dir,'eddy_detections_validation.mat');
    options=struct('make_figure',false,'save_output',true);
    eddy_detect(input_file,detection_file,temporary_dir,options);
    options=struct('make_figure',false,'save_output',false);
    result.eddy_tracking=first_birth(detection_file, ...
        fullfile(temporary_dir,'unused.mat'),temporary_dir,options);
    fprintf('Module 02 validation complete: %d tracks accepted.\n', ...
        nnz(result.eddy_tracking.accepted_track));
else
    result.eddy_tracking=[];
    fprintf('Module 02 skipped; install required functions: %s\n', ...
        strjoin(eddy_missing,', '));
end

velocity_file=fullfile(root,'03_barotropic_instability','data','reanaly_58_movavg.mat');
topography_file=fullfile(root,'03_barotropic_instability','data', ...
    'smoothed_topography_on_xg.mat');
options=struct('wavelength_min_km',10,'wavelength_max_km',300, ...
    'wavelength_step_km',10,'chunk_size',5,'use_parallel',false, ...
    'resume',false,'make_figure',false,'save_output',false);
result.barotropic=time_var_reanalysis_58(velocity_file,topography_file,'',options);
fprintf('Module 03 validation complete: %d profiles processed.\n', ...
    numel(result.barotropic.time));

required={'sw_f','sw_dist','sw_pden','othercolor'};
missing=required(cellfun(@(f) exist(f,'file')~=2,required));
if isempty(missing)
    options=struct('save_output',false,'make_figure',false);
    result.baroclinic=baroclinic_instability_odv_58N22p5W([],options);
    fprintf('Module 04 validation complete at %.1f W, %.1f N.\n', ...
        abs(result.baroclinic.longitude),result.baroclinic.latitude);
else
    result.baroclinic=[];
    fprintf('Module 04 skipped; install source-program functions: %s\n', ...
        strjoin(missing,', '));
end
end

function remove_validation_temp(path_value)
if exist(path_value,'dir')==7, rmdir(path_value,'s'); end
end
