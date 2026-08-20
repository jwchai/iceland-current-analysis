%% Track anticyclonic eddies and calculate their accepted birth locations
% Run eddy_detect.m first. Paths are relative to the package root.
clearvars; close all; clc;

script_dir=fileparts(mfilename('fullpath'));
package_root=fileparts(script_dir);
oceaneddies_dir=fullfile(package_root,'external','OceanEddies');
detection_dir=fullfile(package_root,'outputs','eddy_detection');
work_dir=fullfile(package_root,'outputs','eddy_tracking');
split_dir=fullfile(work_dir,'split','anticyclonic');

type='anticyclonic';
minimum_surface_area_m2=8e9;
time_frequency_days=1;
tracking_tolerance_frames=3;
tracking_minimum_pixels=6;
minimum_track_days=27;
minimum_residence_days=8;
lon_polygon=[-33.8 -23.5 -13.8 -23.8 -33.8];
lat_polygon=[58.2 63.2 60.0 54.0 58.2];
overwrite=false;

if exist(detection_dir,'dir')~=7, error('Missing detection output: %s',detection_dir); end
if exist(oceaneddies_dir,'dir')~=7, error('Missing OceanEddies: %s',oceaneddies_dir); end
addpath(genpath(oceaneddies_dir));
if exist('tolerance_track_lnn','file')~=2, error('tolerance_track_lnn.m is unavailable.'); end
if exist(split_dir,'dir')~=7, mkdir(split_dir); end

if exist(split_dir,'dir')==7
    existing=dir(split_dir); existing=existing(~ismember({existing.name},{'.','..'}));
    if ~isempty(existing)&&~overwrite
        error('Split directory is not empty: %s. Use a clean directory or set overwrite=true.',split_dir);
    elseif ~isempty(existing)
        rmdir(split_dir,'s');
        mkdir(split_dir);
    end
end

%% Split the combined detection files and retain the original area criterion
files=dir(fullfile(detection_dir,'eddies_*.mat'));
if isempty(files), error('No eddies_*.mat files found in %s.',detection_dir); end
dates=NaT(numel(files),1); source=cell(numel(files),1);
for k=1:numel(files)
    S=load(fullfile(files(k).folder,files(k).name),'eddies'); E=S.eddies;
    if ~isfield(E,'time_datetime')||~isdatetime(E.time_datetime)
        error('Detection file lacks corrected time_datetime: %s',files(k).name);
    end
    dates(k)=E.time_datetime; source{k}=fullfile(files(k).folder,files(k).name);
end
[dates,order]=sort(dates); source=source(order);
if numel(unique(dates))~=numel(dates)||any(days(diff(dates))~=time_frequency_days)
    error('Detection files must contain one unique file for every consecutive day.');
end

% The upstream tracker cannot load an entirely empty daily file. Therefore,
% days with no eddy passing the physical-area threshold are treated as hard
% gaps: the remaining dates are separated into continuous tracking segments.
segment_id=0; previous_kept_index=NaN;
segment_first_index=[]; segment_last_index=[];
for k=1:numel(source)
    S=load(source{k},'eddies'); candidates=S.eddies.anticyc;
    if isempty(candidates), continue; end
    if ~isfield(candidates,'SurfaceArea')
        error('SurfaceArea is missing on %s.',string(dates(k),'yyyy-MM-dd'));
    end
    eddies=candidates([candidates.SurfaceArea]>=minimum_surface_area_m2);
    if isempty(eddies), continue; end
    if isnan(previous_kept_index)||k~=previous_kept_index+1
        segment_id=segment_id+1;
        segment_first_index(segment_id)=k; %#ok<SAGROW>
    end
    segment_last_index(segment_id)=k; %#ok<SAGROW>
    segment_dir=fullfile(split_dir,sprintf('segment_%04d',segment_id));
    if exist(segment_dir,'dir')~=7, mkdir(segment_dir); end
    save(fullfile(segment_dir,sprintf('anticyclonic_%s.mat', ...
        string(dates(k),'yyyyMMdd'))),'eddies');
    previous_kept_index=k;
end
if segment_id==0, error('No anticyclonic eddy passes the surface-area threshold.'); end

%% OceanEddies tolerance tracking, independently for each continuous segment
tracks={}; revived=0; dropped_criteria=0; dropped_vicinity=0;
for s=1:segment_id
    if segment_last_index(s)-segment_first_index(s)+1<2, continue; end
    segment_dir=fullfile(split_dir,sprintf('segment_%04d',s));
    [segment_tracks,r,dc,dv]=tolerance_track_lnn(segment_dir,type, ...
        time_frequency_days,tracking_tolerance_frames,tracking_minimum_pixels);
    % Column 3 is the tracker-local time-step index. Convert it to the
    % corresponding index in the complete, sorted daily date vector.
    for j=1:numel(segment_tracks)
        segment_tracks{j}(:,3)=segment_tracks{j}(:,3)+segment_first_index(s)-1;
    end
    tracks=[tracks segment_tracks]; %#ok<AGROW>
    revived=revived+r; dropped_criteria=dropped_criteria+dc;
    dropped_vicinity=dropped_vicinity+dv;
end

%% Birth, lifetime, and residence criteria using the real daily dates
birth_time=NaT(numel(tracks),1); birth_lon=nan(numel(tracks),1);
birth_lat=nan(numel(tracks),1); accepted=false(numel(tracks),1);
track_duration_days=nan(numel(tracks),1); residence_days=nan(numel(tracks),1);
n_accepted=0;
for k=1:numel(tracks)
    M=tracks{k};
    if isempty(M)||size(M,1)<2||size(M,2)<3, continue; end
    step=round(M(:,3));
    if any(step<1|step>numel(dates)), error('Track %d contains an invalid time-step index.',k); end
    track_dates=dates(step);
    duration=days(track_dates(end)-track_dates(1))+time_frequency_days;
    track_duration_days(k)=duration;
    if duration<minimum_track_days, continue; end
    lat=double(M(:,1)); lon=double(M(:,2));
    if ~inpolygon(lon(1),lat(1),lon_polygon,lat_polygon), continue; end
    inside=inpolygon(lon,lat,lon_polygon,lat_polygon);
    residence=sum(inside)*time_frequency_days; residence_days(k)=residence;
    if residence<minimum_residence_days, continue; end
    accepted(k)=true;
    n_accepted=n_accepted+1;
    birth_lat(n_accepted)=lat(1); birth_lon(n_accepted)=lon(1);
    birth_time(n_accepted)=track_dates(1);
end
birth_lat=birth_lat(1:n_accepted); birth_lon=birth_lon(1:n_accepted);
birth_time=birth_time(1:n_accepted);

%% Birth probability per grid cell and per physical area
lon_edges=-40:0.2:-12; lat_edges=52:0.2:65;
[counts,lat_edges,lon_edges]=histcounts2(birth_lat,birth_lon,lat_edges,lon_edges);
if sum(counts(:))>0, probability_mass=counts/sum(counts(:));
else, probability_mass=zeros(size(counts));
end
R_km=6371; dlon=deg2rad(diff(lon_edges));
area_km2=zeros(numel(lat_edges)-1,numel(lon_edges)-1);
for j=1:numel(lat_edges)-1
    band=R_km^2*abs(sind(lat_edges(j+1))-sind(lat_edges(j)));
    area_km2(j,:)=band*dlon;
end
probability_density_km2=probability_mass./area_km2;

result=struct('tracks',{tracks},'accepted_track',accepted,'birth_lon',birth_lon, ...
    'birth_lat',birth_lat,'birth_time',birth_time,'track_duration_days', ...
    track_duration_days,'residence_days',residence_days,'counts',counts, ...
    'probability_mass',probability_mass,'probability_density_per_km2', ...
    probability_density_km2,'latitude_edges',lat_edges,'longitude_edges',lon_edges, ...
    'revived',revived,'dropped_criteria',dropped_criteria, ...
    'dropped_vicinity',dropped_vicinity,'minimum_surface_area_m2', ...
    minimum_surface_area_m2,'tracking_tolerance_frames',tracking_tolerance_frames, ...
    'tracking_minimum_pixels',tracking_minimum_pixels);
if exist(work_dir,'dir')~=7, mkdir(work_dir); end
save(fullfile(work_dir,'eddy_tracking_and_births.mat'),'result','-v7.3');

%% Plot normalized birth probability per grid cell
lon_centres=lon_edges(1:end-1)+diff(lon_edges)/2;
lat_centres=lat_edges(1:end-1)+diff(lat_edges)/2;
figure('Color','w');
imagesc(lon_centres,lat_centres,probability_mass); set(gca,'YDir','normal');
hold on; plot(lon_polygon,lat_polygon,'k-','LineWidth',1.5);
xlabel('Longitude'); ylabel('Latitude'); colorbar;
title('Anticyclonic eddy birth probability per grid cell');
exportgraphics(gcf,fullfile(work_dir,'anticyclonic_birth_probability.png'), ...
    'Resolution',300,'BackgroundColor','white');
fprintf('Accepted %d tracks. Results: %s\n',nnz(accepted),work_dir);
