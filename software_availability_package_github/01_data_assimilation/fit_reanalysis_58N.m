function result=fit_reanalysis_58N(data_dir,output_dir,options)
%FIT_REANALYSIS_58N Mooring-constrained GLORYS velocity fit along 58 N.
%
% result = fit_reanalysis_58N(data_dir,output_dir,options)
%
% This is the function version of the author's original
% fit_reanaylse_for_fig2.m workflow. It reads the two GLORYS sections and
% the four UMM2/IB3/IB4/IB5 mooring records used in the manuscript, performs
% the weighted least-squares fit, interpolates the result to the 1-km grid,
% and calculates the 30-day background current used by the instability
% analysis. Both GLORYS and the mooring observations are represented as
% nominal 0--215 m depth averages. A centered 14-day ADCP mean is calculated
% before each daily reconstruction. No post-hoc velocity scaling is applied.
%
% Example
%   data_dir = '/path/to/raw_input_data';
%   output_dir = '/path/to/output';
%   result = fit_reanalysis_58N(data_dir,output_dir,struct());
%
% options fields and manuscript defaults
%   adcp_filter_window_days = 14  (total width; +/-7 days)
%   gamma                   = 400
%   lambda_background      = 0.5
%   lambda_smoothness      = 50
%   background_window_days = 30
%   make_figure            = true
%   figure_index           = 500
%   save_output            = true

if nargin<1 || isempty(data_dir), error('data_dir is required.'); end
if nargin<2 || isempty(output_dir), output_dir=''; end
if nargin<3 || isempty(options), options=struct(); end
opt=defaults(options);
if exist(data_dir,'dir')~=7, error('Data directory not found: %s',data_dir); end
if opt.save_output && isempty(output_dir), error('output_dir is required when save_output=true.'); end
if ~isempty(output_dir) && exist(output_dir,'dir')~=7, mkdir(output_dir); end

%% 1. The two GLORYS daily files used in the original analysis
glorys_files={ ...
    fullfile(data_dir,'cmems_mod_glo_phy_my_0.083deg_P1D-m_1755700424406.nc'), ...
    fullfile(data_dir,'cmems_mod_glo_phy_myint_0.083deg_P1D-m_1756090265159.nc')};
for j=1:numel(glorys_files)
    if exist(glorys_files{j},'file')~=2, error('Missing GLORYS file: %s',glorys_files{j}); end
end

[time1,lon1,v1]=read_glorys_velocity(glorys_files{1},opt.depth_range_m);
[time2,lon2,v2]=read_glorys_velocity(glorys_files{2},opt.depth_range_m);
if numel(lon1)~=numel(lon2) || any(abs(lon1-lon2)>1e-10)
    error('The two GLORYS longitude grids are not identical.');
end
time_bg=[time1;time2]; v_bg=[v1;v2]; lon_bg=lon1(:);
[time_bg,order]=sort(time_bg); v_bg=v_bg(order,:);
[time_bg,unique_index]=unique(time_bg,'stable'); v_bg=v_bg(unique_index,:);

%% 2. Nominal 0--215 m averages from the four 58 N moorings
lon_obs=[-28;-24.42;-24.14;-19.14];
umm2_file1=locate_umm2_file(data_dir,'OS_OSNAP-UMM2_201807_ADCP_300m.nc');
umm2_file2=locate_umm2_file(data_dir,'OS_OSNAP-UMM2_202009_ADCP_300m.nc');
[time_uum2,v_uum2]=load_umm2_raw_depth_average( ...
    umm2_file1,umm2_file2,opt.depth_range_m);
[time_ib3,v_ib3]=load_uk_depth_average( ...
    fullfile(data_dir,'2018-2020 UK ib3 0-100 .mat'), ...
    fullfile(data_dir,'2020-2022 UK ib3 0-100 .mat'), ...
    fullfile(data_dir,'2018-2020 UK ib3 200-100 .mat'), ...
    fullfile(data_dir,'2020-2022 UK ib3 200-100 .mat'),opt.mooring_layer_thickness_m);
[time_ib4,v_ib4]=load_uk_depth_average( ...
    fullfile(data_dir,'2018-2020 UK ib4 0-100 .mat'), ...
    fullfile(data_dir,'2020-2022 UK ib4 0-100 .mat'), ...
    fullfile(data_dir,'2018-2020 UK ib4 200-100 .mat'), ...
    fullfile(data_dir,'2020-2022 UK ib4 200-100 .mat'),opt.mooring_layer_thickness_m);
[time_ib5,v_ib5]=load_uk_depth_average( ...
    fullfile(data_dir,'2018-2020 UK ib5 0-100 .mat'), ...
    fullfile(data_dir,'2020-2022 UK ib5 0-100 .mat'), ...
    fullfile(data_dir,'2018-2020 UK ib5 200-100 .mat'), ...
    fullfile(data_dir,'2020-2022 UK ib5 200-100 .mat'),opt.mooring_layer_thickness_m);
time_cell={time_uum2,time_ib3,time_ib4,time_ib5};
velocity_cell={v_uum2,v_ib3,v_ib4,v_ib5};

start_common=max(cellfun(@(t) min(t,[],'omitnan'),time_cell));
end_common=min(cellfun(@(t) max(t,[],'omitnan'),time_cell));
if start_common>=end_common, error('The four mooring records have no common period.'); end

%% 3. Mooring-constrained weighted least-squares reconstruction
nlon=numel(lon_bg); ntime=numel(time_bg); nmoor=numel(lon_obs);
e=ones(nlon,1); D2=spdiags([e -2*e e],[-1 0 1],nlon,nlon);
D2(1,:)=0; D2(end,:)=0;
H=zeros(nmoor,nlon);
for j=1:nmoor
    if lon_obs(j)<min(lon_bg) || lon_obs(j)>max(lon_bg)
        error('Mooring longitude %.2f is outside the GLORYS section.',lon_obs(j));
    end
    H(j,:)=interp1(lon_bg,eye(nlon),lon_obs(j),'linear',0);
end

lat0=58; lon_start=-28; xg=0:1:520;
lon_xg=lon_start+xg/(111.32*cosd(lat0));
v_fitted=nan(ntime,nlon); v_on_1km=nan(ntime,numel(xg));
observations=nan(ntime,nmoor); fit_residual_rms=nan(ntime,1);
half_window=opt.adcp_filter_window_days/2;
A=[sqrt(opt.gamma)*H; ...
   sqrt(opt.lambda_background)*eye(nlon); ...
   sqrt(opt.lambda_smoothness)*D2];

for it=1:ntime
    for j=1:nmoor
        % This is the centered 14-day ADCP low-pass preprocessing described
        % in the manuscript, evaluated at each daily GLORYS analysis time.
        use=abs(time_cell{j}-time_bg(it))<=half_window;
        values=velocity_cell{j}(use); values=values(isfinite(values));
        if ~isempty(values), observations(it,j)=mean(values); end
    end
    if all(isfinite(observations(it,:))) && all(isfinite(v_bg(it,:)))
        b=[sqrt(opt.gamma)*observations(it,:)'; ...
           sqrt(opt.lambda_background)*v_bg(it,:)'; ...
           zeros(nlon,1)];
        v_fit=A\b;
        v_fitted(it,:)=v_fit';
        v_on_1km(it,:)=interp1(lon_bg,v_fit,lon_xg,'linear',nan);
        fit_residual_rms(it)=sqrt(mean((H*v_fit-observations(it,:)').^2));
    end
end

%% 4. Common period and 30-day quasi-steady background current
use_common=time_bg>=ceil(start_common) & time_bg<=floor(end_common);
time=time_bg(use_common); vg=v_on_1km(use_common,:);
if isempty(time), error('No reconstructed GLORYS days fall in the common mooring period.'); end
dt_day=median(diff(time),'omitnan');
window_count=max(1,round(opt.background_window_days/dt_day));
if mod(window_count,2)==0, window_count=window_count+1; end
vg_smooth=movmean(vg,window_count,1,'omitnan');
time_day=(ceil(min(time)):floor(max(time)))';
vg_day=nan(numel(time_day),size(vg,2));
for j=1:numel(time_day)
    [~,nearest]=min(abs(time-time_day(j))); vg_day(j,:)=vg_smooth(nearest,:);
end

result=struct('latitude',lat0,'longitude',lon_bg,'background_velocity',v_bg, ...
    'background_time',time_bg,'mooring_longitude',lon_obs, ...
    'mooring_time',{time_cell},'mooring_velocity',{velocity_cell}, ...
    'mooring_windowed_velocity',observations,'fitted_velocity',v_fitted, ...
    'fit_residual_rms',fit_residual_rms,'xg_km',xg,'lon_xg',lon_xg, ...
    'time',time,'vg',vg,'vg_smooth',vg_smooth,'time_day',time_day, ...
    'vg_day',vg_day,'depth_range_m',opt.depth_range_m, ...
    'adcp_filter_window_days',opt.adcp_filter_window_days, ...
    'options',opt,'source_files',{glorys_files});

%% 5. Save the variables expected by the downstream 58 N analysis
if opt.save_output
    save(fullfile(output_dir,'reanaly_fit_58N.mat'),'xg','lon_xg','vg','time', ...
        'vg_smooth','time_day','vg_day','window_count','opt','-v7.3');
end

if opt.make_figure
    available=find(all(isfinite(v_fitted),2));
    if isempty(available), error('No complete fitted profile is available for plotting.'); end
    requested=min(max(1,round(opt.figure_index)),ntime);
    [~,q]=min(abs(available-requested)); iday=available(q);
    figure('Color','w');
    plot(lon_bg,v_bg(iday,:),'k--','LineWidth',1.2); hold on;
    plot(lon_bg,v_fitted(iday,:),'b-','LineWidth',2);
    plot(lon_obs,observations(iday,:),'ro','MarkerSize',7,'LineWidth',1.3);
    xlabel('Longitude (degree E)'); ylabel('Meridional velocity (m s^{-1})');
    legend('GLORYS background','Mooring-constrained fit','Mooring observations','Location','best');
    fit_date=datetime(time_bg(iday),'ConvertFrom','datenum','Format','yyyy-MM-dd');
    title(sprintf('58 N fit: %s',char(fit_date))); grid on;
    if opt.save_output
        exportgraphics(gcf,fullfile(output_dir,'reanalysis_fit_58N.png'),'Resolution',300);
    end
end
end

function opt=defaults(opt)
d=struct('adcp_filter_window_days',14,'gamma',400,'lambda_background',0.5, ...
    'lambda_smoothness',50,'background_window_days',30,'make_figure',true, ...
    'figure_index',500,'save_output',true,'depth_range_m',[0 215], ...
    'mooring_layer_thickness_m',[100 115]);
names=fieldnames(d);
for j=1:numel(names)
    if ~isfield(opt,names{j}) || isempty(opt.(names{j})), opt.(names{j})=d.(names{j}); end
end
if ~isscalar(opt.adcp_filter_window_days) || opt.adcp_filter_window_days<=0
    error('adcp_filter_window_days must be a positive scalar.');
end
if ~isequal(size(opt.depth_range_m),[1 2]) || opt.depth_range_m(1)~=0 || ...
        opt.depth_range_m(2)~=215
    error('The manuscript configuration requires depth_range_m=[0 215].');
end
if numel(opt.mooring_layer_thickness_m)~=2 || ...
        any(opt.mooring_layer_thickness_m<=0) || ...
        abs(sum(opt.mooring_layer_thickness_m)-215)>1e-10
    error('mooring_layer_thickness_m must contain two positive values summing to 215 m.');
end
end

function [time,lon,vmean]=read_glorys_velocity(file,depth_range_m)
lon=double(ncread(file,'longitude')); lon=lon(:);
depth=double(ncread(file,'depth')); depth=depth(:);
raw_time=double(ncread(file,'time')); raw_time=raw_time(:);
units=ncreadatt(file,'time','units'); time=cf_time_to_datenum(raw_time,units);
raw=squeeze(double(ncread(file,'vo')));
sz=size(raw); target=[numel(time),numel(depth),numel(lon)];
perm=zeros(1,3); remaining=1:3;
for j=1:3
    candidate=remaining(sz(remaining)==target(j));
    if isempty(candidate), error('Cannot identify vo dimensions in %s.',file); end
    perm(j)=candidate(1); remaining(remaining==perm(j))=[];
end
raw=permute(raw,perm);
% Finite-volume layer-overlap weights give a true thickness-weighted
% average over 0--215 m, rather than an arithmetic mean over all 50 levels.
edges=[0;0.5*(depth(1:end-1)+depth(2:end));Inf];
weights=max(0,min(edges(2:end),depth_range_m(2))-max(edges(1:end-1),depth_range_m(1)));
if abs(sum(weights)-diff(depth_range_m))>1e-8
    error('GLORYS depth coordinate does not cover the requested 0--215 m layer.');
end
w=reshape(weights,1,[],1);
valid=isfinite(raw);
numerator=sum(raw.*w,2,'omitnan');
denominator=sum(valid.*w,2);
vmean=squeeze(numerator./denominator);
vmean(squeeze(denominator)<0.95*diff(depth_range_m))=NaN;
if ~isequal(size(vmean),[numel(time),numel(lon)]), error('Unexpected depth-mean vo dimensions.'); end
end

function dn=cf_time_to_datenum(value,units)
tok=regexp(strtrim(units),'^(seconds|hours|days)\s+since\s+(.+)$','tokens','once','ignorecase');
if isempty(tok), error('Unsupported NetCDF time units: %s',units); end
origin_text=regexprep(strrep(strtrim(tok{2}),'T',' '),'\s*(UTC|Z)$','');
origin=datenum(origin_text); %#ok<DATNM>
switch lower(tok{1})
    case 'seconds', dn=origin+value/86400;
    case 'hours', dn=origin+value/24;
    case 'days', dn=origin+value;
end
end

function [time,velocity]=load_mooring_pair(file1,file2,velocity_name)
if exist(file1,'file')~=2 || exist(file2,'file')~=2, error('Missing mooring file pair.'); end
S1=load(file1); S2=load(file2);
if ~isfield(S1,'time') || ~isfield(S2,'time') || ...
        ~isfield(S1,velocity_name) || ~isfield(S2,velocity_name)
    error('Expected time and %s in mooring files.',velocity_name);
end
time=[double(S1.time(:));double(S2.time(:))];
v1=S1.(velocity_name); v2=S2.(velocity_name);
velocity=[double(v1(:));double(v2(:))];
valid=isfinite(time); time=time(valid); velocity=velocity(valid);
[time,order]=sort(time); velocity=velocity(order);
[time,unique_index]=unique(time,'stable'); velocity=velocity(unique_index);
end

function file=locate_umm2_file(data_dir,name)
candidates={fullfile(data_dir,name),fullfile(data_dir,'UMM2',name), ...
    fullfile(fileparts(data_dir),'UMM2',name)};
found=find(cellfun(@(p) exist(p,'file')==2,candidates),1);
if isempty(found)
    error('Missing raw UMM2 NetCDF %s. Place it in data_dir or data_dir/UMM2.',name);
end
file=candidates{found};
end

function [time,velocity]=load_umm2_raw_depth_average(file1,file2,depth_range_m)
[t1,v1]=load_umm2_raw_single(file1,depth_range_m);
[t2,v2]=load_umm2_raw_single(file2,depth_range_m);
time=[t1;t2]; velocity=[v1;v2];
[time,order]=sort(time); velocity=velocity(order);
[time,unique_index]=unique(time,'stable'); velocity=velocity(unique_index);
end

function [time,velocity]=load_umm2_raw_single(file,depth_range_m)
if exist(file,'file')~=2, error('Missing raw UMM2 file: %s',file); end
raw_time=double(ncread(file,'TIME')); raw_time=raw_time(:);
units=ncreadatt(file,'TIME','units'); time=cf_time_to_datenum(raw_time,units);
depth=double(ncread(file,'BINDEPTH'));
v=double(ncread(file,'VCUR'));
if size(v,2)==numel(time)
    % MATLAB commonly returns these NetCDF variables as [bin x time].
elseif size(v,1)==numel(time)
    v=v'; depth=depth';
else
    error('Cannot identify TIME dimension in %s.',file);
end
if ~isequal(size(v),size(depth)), error('VCUR and BINDEPTH sizes differ in %s.',file); end
v(abs(v)>10)=NaN; depth(depth<0 | depth>12000)=NaN;
velocity=nan(numel(time),1);
for it=1:numel(time)
    use=depth(:,it)>=depth_range_m(1) & depth(:,it)<=depth_range_m(2) & ...
        isfinite(v(:,it));
    if nnz(use)>=2, velocity(it)=mean(v(use,it),'omitnan'); end
end
end

function [time,velocity]=load_uk_depth_average(upper1,upper2,lower1,lower2,layer_thickness)
[time_upper,velocity_upper]=load_mooring_pair(upper1,upper2,'v_mean');
[time_lower,velocity_lower]=load_mooring_pair(lower1,lower2,'v_mean');
time=time_upper;
if numel(time_lower)<2, error('Too few lower-layer samples in UK mooring files.'); end
nearest=interp1(time_lower,(1:numel(time_lower))',time,'nearest',NaN);
matched=nan(size(time)); usable=isfinite(nearest);
nearest_index=round(nearest(usable));
dt=max(1e-6,0.51*median(diff(time_lower),'omitnan'));
usable_index=find(usable);
close_enough=abs(time_lower(nearest_index)-time(usable))<=dt;
usable(:)=false; usable(usable_index(close_enough))=true;
nearest_index=round(nearest(usable));
matched(usable)=velocity_lower(nearest_index);
velocity=nan(size(time)); complete=isfinite(velocity_upper)&isfinite(matched);
velocity(complete)=(layer_thickness(1)*velocity_upper(complete)+ ...
    layer_thickness(2)*matched(complete))/sum(layer_thickness);
end
