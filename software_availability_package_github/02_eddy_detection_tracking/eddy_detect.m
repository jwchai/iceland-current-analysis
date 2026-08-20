%% Detect cyclonic and anticyclonic eddies from the CMEMS ADT product
% Run this script from any working directory. All paths are resolved from
% the software-availability package root.
clearvars; close all; clc;

script_dir=fileparts(mfilename('fullpath'));
package_root=fileparts(script_dir);
data_dir=fullfile(package_root,'data');
oceaneddies_dir=fullfile(package_root,'external','OceanEddies');
output_dir=fullfile(package_root,'outputs','eddy_detection');

input_name= ...
    'cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.125deg_P1D_adt-sla_39.94W-0.06W_50.06N-64.94N_1993-01-01-2025-10-18.nc';
input_file=fullfile(data_dir,input_name);
ssh_variable='adt'; longitude_variable='longitude';
latitude_variable='latitude'; time_variable='time';
scan_version='v2'; detection_minimum_pixels=9;
is_padding=false; % regional domain, not a global periodic longitude grid
overwrite=false;

if exist(input_file,'file')~=2, error('Missing CMEMS input: %s',input_file); end
if exist(oceaneddies_dir,'dir')~=7
    error('Missing OceanEddies dependency: %s',oceaneddies_dir);
end
addpath(genpath(oceaneddies_dir));
if exist('scan_single','file')~=2, error('scan_single.m is not on the MATLAB path.'); end
if exist(output_dir,'dir')~=7, mkdir(output_dir); end

lat=double(ncread(input_file,latitude_variable)); lat=lat(:);
lon=double(ncread(input_file,longitude_variable)); lon=lon(:);
lon(lon>=180)=lon(lon>=180)-360;
raw_time=double(ncread(input_file,time_variable)); raw_time=raw_time(:);
time_units=ncreadatt(input_file,time_variable,'units');
try
    calendar=ncreadatt(input_file,time_variable,'calendar');
catch
    calendar='standard';
end
time_datetime=decode_cf_time(raw_time,time_units,calendar);

area_map=grid_cell_area(lat,lon);
vinfo=ncinfo(input_file,ssh_variable);
dim_names=string({vinfo.Dimensions.Name});
lon_dim=find(dim_names==longitude_variable,1);
lat_dim=find(dim_names==latitude_variable,1);
time_dim=find(dim_names==time_variable,1);
if isempty(lon_dim)||isempty(lat_dim)||isempty(time_dim)
    error('Cannot identify longitude, latitude, and time dimensions of %s.',ssh_variable);
end

for it=1:numel(raw_time)
    start=ones(1,numel(vinfo.Size)); count=vinfo.Size;
    start(time_dim)=it; count(time_dim)=1;
    ssh=squeeze(double(ncread(input_file,ssh_variable,start,count)));
    if isequal(size(ssh),[numel(lon),numel(lat)])
        ssh=ssh';
    elseif ~isequal(size(ssh),[numel(lat),numel(lon)])
        error('Unexpected ADT slice size at time index %d.',it);
    end
    ssh(~isfinite(ssh))=NaN;
    dn=datenum(time_datetime(it)); %#ok<DATNM>
    cyclonic=scan_single(ssh,lat,lon,dn,'cyclonic',scan_version,area_map, ...
        'sshUnits','meters','minimumArea',detection_minimum_pixels, ...
        'isPadding',is_padding);
    anticyc=scan_single(ssh,lat,lon,dn,'anticyc',scan_version,area_map, ...
        'sshUnits','meters','minimumArea',detection_minimum_pixels, ...
        'isPadding',is_padding);
    eddies=struct('cyclonic',cyclonic,'anticyc',anticyc,'lat',lat,'lon',lon, ...
        'time_datetime',time_datetime(it),'time_datenum',dn, ...
        'time_raw',raw_time(it),'time_units',time_units,'calendar',calendar, ...
        'source_file_relative',fullfile('data',input_name),'scan_version',scan_version, ...
        'detection_minimum_pixels',detection_minimum_pixels,'is_padding',is_padding);
    output_file=fullfile(output_dir,sprintf('eddies_%s.mat', ...
        datestr(time_datetime(it),'yyyymmdd'))); %#ok<DATST>
    if exist(output_file,'file')==2 && ~overwrite
        error('Output already exists: %s. Set overwrite=true to replace it.',output_file);
    end
    save(output_file,'eddies','-v7.3');
    if mod(it,20)==0||it==numel(raw_time), fprintf('Saved %d / %d\n',it,numel(raw_time)); end
end
fprintf('Eddy detection completed: %s\n',output_dir);

%% Optional view of the final time slice
figure('Color','w');
h1=imagesc(lon,lat,ssh); set(gca,'YDir','normal'); hold on;
if isempty(anticyc), h2=scatter(nan,nan,12,'r','filled');
else, h2=scatter([anticyc.Lon],[anticyc.Lat],12,'r','filled'); end
if isempty(cyclonic), h3=scatter(nan,nan,12,'b','filled');
else, h3=scatter([cyclonic.Lon],[cyclonic.Lat],12,'b','filled'); end
xlabel('Longitude'); ylabel('Latitude'); colorbar;
legend([h1 h2 h3],{'ADT','Anticyclonic centre','Cyclonic centre'},'Location','best');
date_label=string(time_datetime(end),'yyyy-MM-dd');
title(sprintf('OceanEddies detection, %s',date_label));

function dt=decode_cf_time(value,units,calendar)
if ~any(strcmpi(calendar,{'standard','gregorian','proleptic_gregorian'}))
    error('Unsupported NetCDF calendar: %s',calendar);
end
token=regexp(strtrim(units),'^(seconds|hours|days)\s+since\s+(.+)$', ...
    'tokens','once','ignorecase');
if isempty(token), error('Unsupported NetCDF time units: %s',units); end
origin_text=regexprep(strrep(strtrim(token{2}),'T',' '),'\s*(UTC|Z)$','');
origin=datetime(origin_text,'TimeZone','UTC');
switch lower(token{1})
    case 'seconds', dt=origin+seconds(value);
    case 'hours', dt=origin+hours(value);
    case 'days', dt=origin+days(value);
end
dt.TimeZone='';
end

function area=grid_cell_area(lat,lon)
R=6371000;
dlat=abs(median(diff(lat))); dlon=abs(median(diff(lon)))*pi/180;
area=zeros(numel(lat),numel(lon));
for j=1:numel(lat)
    south=(lat(j)-dlat/2)*pi/180; north=(lat(j)+dlat/2)*pi/180;
    area(j,:)=R^2*dlon*abs(sin(north)-sin(south));
end
end
