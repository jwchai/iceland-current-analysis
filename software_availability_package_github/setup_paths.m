function repository_root=setup_paths()
%SETUP_PATHS Add repository code and installed dependencies to MATLAB path.

repository_root=fileparts(mfilename('fullpath'));
module_dirs={ ...
    fullfile(repository_root,'01_data_assimilation'), ...
    fullfile(repository_root,'02_eddy_detection_tracking'), ...
    fullfile(repository_root,'03_barotropic_instability'), ...
    fullfile(repository_root,'04_baroclinic_instability')};
for j=1:numel(module_dirs)
    if exist(module_dirs{j},'dir')==7, addpath(module_dirs{j}); end
end

dependency_dirs={ ...
    fullfile(repository_root,'external','OceanEddies'), ...
    fullfile(repository_root,'external','seawater')};
for j=1:numel(dependency_dirs)
    if exist(dependency_dirs{j},'dir')==7
        addpath(genpath(dependency_dirs{j}));
    end
end
end

