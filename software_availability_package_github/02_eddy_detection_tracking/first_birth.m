function result=first_birth(detection_file,output_file,figure_dir,options)
%FIRST_BIRTH Track compact-example anticyclones and map their births.
%
% The default thresholds are deliberately relaxed so that the bundled
% 15-day sample validates the complete detection-to-tracking workflow. See
% README.md for the stricter manuscript-analysis thresholds.

module_dir=fileparts(mfilename('fullpath'));
if nargin<1||isempty(detection_file)
    detection_file=fullfile(module_dir,'data','eddy_detections_validation.mat');
end
if nargin<2||isempty(output_file)
    output_file=fullfile(module_dir,'data','eddy_tracking_validation.mat');
end
if nargin<3||isempty(figure_dir)
    figure_dir=fullfile(module_dir,'figure');
end
if nargin<4||isempty(options), options=struct(); end
options=defaults(options,struct( ...
    'minimum_surface_area_m2',0, ...
    'time_frequency_days',1, ...
    'tracking_tolerance_frames',1, ...
    'tracking_minimum_pixels',1, ...
    'minimum_track_days',4, ...
    'minimum_residence_days',1, ...
    'longitude_polygon',[-33.8 -23.5 -13.8 -23.8 -33.8], ...
    'latitude_polygon',[58.2 63.2 60.0 54.0 58.2], ...
    'make_figure',true, ...
    'save_output',true));

if exist(detection_file,'file')~=2
    error('Missing detection result: %s. Run eddy_detect first.',detection_file);
end
tracker_file=which('tolerance_track_lnn');
if isempty(tracker_file)
    error(['OceanEddies tolerance_track_lnn.m is unavailable. Install the ' ...
        'dependency under ../external/OceanEddies and run setup_paths.']);
end
S=load(detection_file,'frames','metadata');
frames=S.frames(:);
dates=vertcat(frames.time_datetime);
if numel(unique(dates))~=numel(dates)|| ...
        any(days(diff(dates))~=options.time_frequency_days)
    error('Detection frames must represent unique consecutive days.');
end

temporary_root=tempname;
mkdir(temporary_root);
cleanup_temp=onCleanup(@() remove_temporary(temporary_root)); %#ok<NASGU>

% The upstream tracker does not accept a completely empty daily file.
% Split nonempty days into consecutive segments in a temporary directory.
segment_id=0; previous_kept_index=NaN;
segment_first_index=[]; segment_last_index=[];
for k=1:numel(frames)
    candidates=frames(k).anticyclonic;
    if isempty(candidates), continue; end
    if ~isfield(candidates,'SurfaceArea')
        error('SurfaceArea is missing on %s.',datestr(dates(k),'yyyy-mm-dd')); %#ok<DATST>
    end
    eddies=candidates([candidates.SurfaceArea]>=options.minimum_surface_area_m2); %#ok<NASGU>
    if isempty(eddies), continue; end
    if isnan(previous_kept_index)||k~=previous_kept_index+1
        segment_id=segment_id+1;
        segment_first_index(segment_id)=k; %#ok<AGROW>
    end
    segment_last_index(segment_id)=k; %#ok<AGROW>
    segment_dir=fullfile(temporary_root,sprintf('segment_%04d',segment_id));
    if exist(segment_dir,'dir')~=7, mkdir(segment_dir); end
    save(fullfile(segment_dir,sprintf('anticyclonic_%s.mat', ...
        datestr(dates(k),'yyyymmdd'))),'eddies','-v7'); %#ok<DATST>
    previous_kept_index=k;
end
if segment_id==0
    error('No anticyclonic eddy passes the lightweight area threshold.');
end

tracks={}; revived=0; dropped_criteria=0; dropped_vicinity=0;
for s=1:segment_id
    if segment_last_index(s)-segment_first_index(s)+1<2, continue; end
    segment_dir=fullfile(temporary_root,sprintf('segment_%04d',s));
    [segment_tracks,r,dc,dv]=portable_track(tracker_file,segment_dir,options);
    for j=1:numel(segment_tracks)
        segment_tracks{j}(:,3)=segment_tracks{j}(:,3)+segment_first_index(s)-1;
    end
    tracks=[tracks segment_tracks]; %#ok<AGROW>
    revived=revived+r; dropped_criteria=dropped_criteria+dc;
    dropped_vicinity=dropped_vicinity+dv;
end

accepted=false(numel(tracks),1);
birth_time=NaT(numel(tracks),1); birth_lon=nan(numel(tracks),1);
birth_lat=nan(numel(tracks),1); track_duration_days=nan(numel(tracks),1);
residence_days=nan(numel(tracks),1); n_accepted=0;
for k=1:numel(tracks)
    M=tracks{k};
    if isempty(M)||size(M,1)<2||size(M,2)<3, continue; end
    step=round(M(:,3));
    if any(step<1|step>numel(dates))
        error('Track %d contains an invalid time-step index.',k);
    end
    track_dates=dates(step);
    duration=days(track_dates(end)-track_dates(1))+options.time_frequency_days;
    track_duration_days(k)=duration;
    if duration<options.minimum_track_days, continue; end
    lat=double(M(:,1)); lon=double(M(:,2));
    if ~inpolygon(lon(1),lat(1),options.longitude_polygon,options.latitude_polygon)
        continue;
    end
    inside=inpolygon(lon,lat,options.longitude_polygon,options.latitude_polygon);
    residence=sum(inside)*options.time_frequency_days;
    residence_days(k)=residence;
    if residence<options.minimum_residence_days, continue; end
    accepted(k)=true; n_accepted=n_accepted+1;
    birth_lat(n_accepted)=lat(1); birth_lon(n_accepted)=lon(1);
    birth_time(n_accepted)=track_dates(1);
end
birth_lat=birth_lat(1:n_accepted); birth_lon=birth_lon(1:n_accepted);
birth_time=birth_time(1:n_accepted);

lon_edges=-40:0.5:-12; lat_edges=52:0.5:65;
[counts,lat_edges,lon_edges]=histcounts2(birth_lat,birth_lon,lat_edges,lon_edges);
if sum(counts(:))>0, probability_mass=counts/sum(counts(:));
else, probability_mass=zeros(size(counts));
end

parameters=options;
result=struct('tracks',{tracks},'accepted_track',accepted, ...
    'birth_lon',birth_lon,'birth_lat',birth_lat,'birth_time',birth_time, ...
    'track_duration_days',track_duration_days,'residence_days',residence_days, ...
    'counts',counts,'probability_mass',probability_mass, ...
    'latitude_edges',lat_edges,'longitude_edges',lon_edges, ...
    'revived',revived,'dropped_criteria',dropped_criteria, ...
    'dropped_vicinity',dropped_vicinity,'parameters',parameters, ...
    'input_metadata',S.metadata, ...
    'sample_only',true, ...
    'sample_note','Relaxed 15-day workflow validation; not manuscript statistics.');

if options.save_output
    ensure_directory(fileparts(output_file));
    save(output_file,'result','-v7');
    fprintf('Tracking result: %s\n',output_file);
end
if options.make_figure
    ensure_directory(figure_dir);
    plot_tracking(result,options,fullfile(figure_dir,'eddy_tracking_validation.png'));
end
fprintf('Tracked %d candidates; accepted %d lightweight-example tracks.\n', ...
    numel(tracks),nnz(accepted));
end

function [tracks,r,dc,dv]=portable_track(tracker_file,segment_dir,options)
original_dir=pwd;
cleanup=onCleanup(@() cd(original_dir)); %#ok<NASGU>
tracker_dir=fileparts(tracker_file);
addpath(tracker_dir);
cd(tracker_dir);
[tracks,r,dc,dv]=tolerance_track_lnn(segment_dir,'anticyclonic', ...
    options.time_frequency_days,options.tracking_tolerance_frames, ...
    options.tracking_minimum_pixels);
end

function plot_tracking(result,options,figure_file)
fig=figure('Color','w','Visible','off','Position',[100 100 1100 500]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile; hold on;
for k=1:numel(result.tracks)
    M=result.tracks{k};
    if isempty(M), continue; end
    if result.accepted_track(k), colour=[0.85 0.1 0.1]; width=1.8;
    else, colour=[0.72 0.72 0.72]; width=0.7; end
    plot(M(:,2),M(:,1),'-','Color',colour,'LineWidth',width);
end
plot(options.longitude_polygon,options.latitude_polygon,'k-','LineWidth',1.3);
if ~isempty(result.birth_lon)
    scatter(result.birth_lon,result.birth_lat,36,'k','filled');
end
xlim([-40 -12]); ylim([52 65]); grid on;
xlabel('Longitude'); ylabel('Latitude');
title(sprintf('Tracks (accepted %d of %d)',nnz(result.accepted_track),numel(result.tracks)));

nexttile;
lon_centres=result.longitude_edges(1:end-1)+diff(result.longitude_edges)/2;
lat_centres=result.latitude_edges(1:end-1)+diff(result.latitude_edges)/2;
imagesc(lon_centres,lat_centres,result.probability_mass); set(gca,'YDir','normal');
hold on; plot(options.longitude_polygon,options.latitude_polygon,'k-','LineWidth',1.3);
xlabel('Longitude'); ylabel('Latitude'); colorbar;
title('Normalized accepted births per grid cell');
exportgraphics(fig,figure_file,'Resolution',200,'BackgroundColor','white');
close(fig);
end

function remove_temporary(path_value)
if exist(path_value,'dir')==7, rmdir(path_value,'s'); end
end

function ensure_directory(path_value)
if exist(path_value,'dir')~=7, mkdir(path_value); end
end

function value=defaults(value,default_value)
names=fieldnames(default_value);
for k=1:numel(names)
    if ~isfield(value,names{k})||isempty(value.(names{k}))
        value.(names{k})=default_value.(names{k});
    end
end
end
