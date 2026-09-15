%% Lightweight module-02 validation: eddy detection, tracking, and figures
% Run this script from any working directory. All paths are module-relative.

module_dir=fileparts(mfilename('fullpath'));
package_root=fileparts(module_dir);
addpath(package_root);
setup_paths();

input_file=fullfile(module_dir,'data','adt_validation_19930101_19930115.nc');
detection_file=fullfile(module_dir,'data','eddy_detections_validation.mat');
tracking_file=fullfile(module_dir,'data','eddy_tracking_validation.mat');
figure_dir=fullfile(module_dir,'figure');

detection=eddy_detect(input_file,detection_file,figure_dir); %#ok<NASGU>
tracking=first_birth(detection_file,tracking_file,figure_dir); %#ok<NASGU>

fprintf('\nModule 02 lightweight validation complete.\n');
fprintf('Data products: %s\n',fullfile(module_dir,'data'));
fprintf('Figures: %s\n',figure_dir);
