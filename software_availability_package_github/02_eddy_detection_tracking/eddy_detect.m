function result=eddy_detect(input_file,output_file,figure_dir,options)
%EDDY_DETECT Detect eddies in the compact CMEMS ADT validation block.
%
% result=eddy_detect() reads data/adt_validation_19930101_19930115.nc,
% runs OceanEddies Eddyscan v2, writes one compact detection file under
% data/, and exports a diagnostic figure under figure/.
%
% This bundled 15-day block is for code-path validation only. It is not the
% full 1993--2025 dataset used for the manuscript statistics.

module_dir=fileparts(mfilename('fullpath'));
if nargin<1||isempty(input_file)
    input_file=fullfile(module_dir,'data','adt_validation_19930101_19930115.nc');
end
if nargin<2||isempty(output_file)
    output_file=fullfile(module_dir,'data','eddy_detections_validation.mat');
end
if nargin<3||isempty(figure_dir)
    figure_dir=fullfile(module_dir,'figure');
end
if nargin<4||isempty(options), options=struct(); end
options=defaults(options,struct( ...
    'ssh_variable','adt', ...
    'longitude_variable','longitude', ...
    'latitude_variable','latitude', ...
    'time_variable','time', ...
    'scan_version','v2', ...
    'detection_minimum_pixels',9, ...
    'is_padding',false, ...
    'make_figure',true, ...
    'save_output',true));

if exist(input_file,'file')~=2
    error('Missing compact ADT input: %s',input_file);
end
scan_file=which('scan_single');
if isempty(scan_file)
    error(['OceanEddies scan_single.m is unavailable. Install the dependency ' ...
        'under ../external/OceanEddies and run setup_paths.']);
end

lat=double(ncread(input_file,options.latitude_variable)); lat=lat(:);
lon=double(ncread(input_file,options.longitude_variable)); lon=lon(:);
lon(lon>=180)=lon(lon>=180)-360;
raw_time=double(ncread(input_file,options.time_variable)); raw_time=raw_time(:);
time_units=ncreadatt(input_file,options.time_variable,'units');
try
    calendar=ncreadatt(input_file,options.time_variable,'calendar');
catch
    calendar='standard';
end
time_datetime=decode_cf_time(raw_time,time_units,calendar);

area_map=grid_cell_area(lat,lon);
vinfo=ncinfo(input_file,options.ssh_variable);
dim_names=string({vinfo.Dimensions.Name});
lon_dim=find(dim_names==options.longitude_variable,1);
lat_dim=find(dim_names==options.latitude_variable,1);
time_dim=find(dim_names==options.time_variable,1);
if isempty(lon_dim)||isempty(lat_dim)||isempty(time_dim)
    error('Cannot identify longitude, latitude, and time dimensions of %s.', ...
        options.ssh_variable);
end

frames=repmat(struct('cyclonic',[],'anticyclonic',[],'time_datetime',NaT, ...
    'time_datenum',NaN),numel(raw_time),1);
last_ssh=[];
for it=1:numel(raw_time)
    start=ones(1,numel(vinfo.Size)); count=vinfo.Size;
    start(time_dim)=it; count(time_dim)=1;
    ssh=squeeze(double(ncread(input_file,options.ssh_variable,start,count)));
    if isequal(size(ssh),[numel(lon),numel(lat)])
        ssh=ssh';
    elseif ~isequal(size(ssh),[numel(lat),numel(lon)])
        error('Unexpected ADT slice size at time index %d.',it);
    end
    ssh(~isfinite(ssh))=NaN;
    dn=datenum(time_datetime(it)); %#ok<DATNM>

    cyclonic=portable_scan(scan_file,ssh,lat,lon,dn,'cyclonic', ...
        options.scan_version,area_map,options);
    anticyclonic=portable_scan(scan_file,ssh,lat,lon,dn,'anticyc', ...
        options.scan_version,area_map,options);
    frames(it).cyclonic=cyclonic;
    frames(it).anticyclonic=anticyclonic;
    frames(it).time_datetime=time_datetime(it);
    frames(it).time_datenum=dn;
    last_ssh=ssh;
    fprintf('Detected eddies for %s (%d/%d).\n', ...
        datestr(time_datetime(it),'yyyy-mm-dd'),it,numel(raw_time)); %#ok<DATST>
end

metadata=struct( ...
    'sample_only',true, ...
    'sample_note','Compact 15-day code-validation subset; not for manuscript statistics.', ...
    'source_file_relative',fullfile('data',string_filename(input_file)), ...
    'longitude',lon,'latitude',lat,'raw_time',raw_time, ...
    'time_units',time_units,'calendar',calendar, ...
    'scan_version',options.scan_version, ...
    'detection_minimum_pixels',options.detection_minimum_pixels, ...
    'is_padding',options.is_padding);

if options.save_output
    ensure_directory(fileparts(output_file));
    save(output_file,'frames','metadata','-v7');
    fprintf('Detection result: %s\n',output_file);
end

figure_file='';
if options.make_figure
    ensure_directory(figure_dir);
    figure_file=fullfile(figure_dir,'eddy_detection_validation.png');
    plot_detection(lon,lat,last_ssh,frames(end),figure_file);
end

result=struct('frames',{frames},'metadata',metadata, ...
    'output_file',output_file,'figure_file',figure_file);
end

function eddies=portable_scan(scan_file,ssh,lat,lon,dn,cyc,version,area_map,options)
original_dir=pwd;
cleanup=onCleanup(@() cd(original_dir)); %#ok<NASGU>
scan_dir=fileparts(scan_file);
addpath(scan_dir);
lib_dir=fullfile(scan_dir,'lib');
if exist(lib_dir,'dir')==7, addpath(lib_dir); end
cd(scan_dir);
eddies=scan_single(ssh,lat,lon,dn,cyc,version,area_map, ...
    'sshUnits','meters', ...
    'minimumArea',options.detection_minimum_pixels, ...
    'isPadding',options.is_padding);
end

function plot_detection(lon,lat,ssh,frame,figure_file)
fig=figure('Color','w','Visible','off','Position',[100 100 900 520]);
h0=imagesc(lon,lat,ssh);
set(h0,'AlphaData',isfinite(ssh));
set(gca,'YDir','normal','Color',[0.88 0.88 0.88]); hold on;
if isempty(frame.anticyclonic)
    h1=scatter(nan,nan,24,'r','filled');
else
    h1=scatter([frame.anticyclonic.Lon],[frame.anticyclonic.Lat],24,'r','filled');
end
if isempty(frame.cyclonic)
    h2=scatter(nan,nan,24,'b','filled');
else
    h2=scatter([frame.cyclonic.Lon],[frame.cyclonic.Lat],24,'b','filled');
end
xlabel('Longitude'); ylabel('Latitude');
cb=colorbar; cb.Label.String='ADT (m)';
legend([h1 h2],{'Anticyclonic centre','Cyclonic centre'},'Location','best');
title(sprintf('OceanEddies validation detection, %s', ...
    datestr(frame.time_datetime,'yyyy-mm-dd'))); %#ok<DATST>
exportgraphics(fig,figure_file,'Resolution',200,'BackgroundColor','white');
close(fig);
end

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

function value=string_filename(path_value)
[~,name,ext]=fileparts(path_value);
value=[name ext];
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
