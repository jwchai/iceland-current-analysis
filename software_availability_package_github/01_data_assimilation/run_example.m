% Run the bundled module-01 assimilation example from any working directory.
script_dir=fileparts(mfilename('fullpath'));
repository_root=fileparts(script_dir);
addpath(repository_root);
setup_paths();

options=struct('save_output',true,'make_figure',true);
result=fit_reanalysis_58N([],script_dir,options);

fprintf('Module 01 example complete: %d days assimilated.\n',numel(result.time));
fprintf('Data:   %s\n',fullfile(script_dir,'data','reanaly_fit_58N.mat'));
fprintf('Figure: %s\n',fullfile(script_dir,'figure','reanalysis_fit_58N.png'));
