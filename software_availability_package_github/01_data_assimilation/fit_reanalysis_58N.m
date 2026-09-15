function result=fit_reanalysis_58N(input_file,output_root,options)
%FIT_REANALYSIS_58N Validate mooring-constrained GLORYS assimilation at 58 N.
%
% result = fit_reanalysis_58N(input_file,output_root,options)
%
% The repository ships a compact, preprocessed validation input instead of
% the multi-year raw GLORYS and mooring files. The input contains:
%   time_bg      [time x 1] daily MATLAB datenums
%   lon_bg       [longitude x 1] GLORYS longitude
%   v_bg         [time x longitude] 0--215 m GLORYS meridional velocity
%   lon_obs      [mooring x 1] mooring longitudes
%   observations [time x mooring] preprocessed mooring velocity
%
% The observations already contain the manuscript preprocessing: nominal
% 0--215 m depth averages and a centered 14-day mean. This function validates
% the assimilation itself, interpolation to 1 km, 30-day background smoothing,
% and plotting. It does not validate raw-data decoding or preprocessing.
%
% Quick validation using the bundled input (plots, but does not save):
%   result = fit_reanalysis_58N();
%
% Save the MAT result and comparison figure:
%   module_dir = fileparts(mfilename('fullpath'));
%   result = fit_reanalysis_58N([],module_dir, ...
%       struct('save_output',true,'figure_date','2020-12-15'));
%
% Main options
%   gamma                   = 400
%   lambda_background      = 0.5
%   lambda_smoothness      = 50
%   background_window_days = 30
%   make_figure            = true
%   figure_date            = [] (middle available day)
%   figure_index           = [] (alternative one-based sample index)
%   save_output            = false

script_dir=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(input_file)
    input_file=fullfile(script_dir,'data','data_assimilation_validation.mat');
end
if nargin<2 || isempty(output_root), output_root=script_dir; end
if nargin<3 || isempty(options), options=struct(); end
opt=defaults(options);
if exist(input_file,'file')~=2, error('Validation input not found: %s',input_file); end
data_output_dir=fullfile(output_root,'data');
figure_output_dir=fullfile(output_root,'figure');
if opt.save_output && exist(data_output_dir,'dir')~=7, mkdir(data_output_dir); end
if opt.save_output && opt.make_figure && exist(figure_output_dir,'dir')~=7
    mkdir(figure_output_dir);
end

%% 1. Load the compact, aligned validation block
S=load(input_file);
required={'time_bg','lon_bg','v_bg','lon_obs','observations'};
for j=1:numel(required)
    if ~isfield(S,required{j}), error('Validation input is missing %s.',required{j}); end
end
time_bg=double(S.time_bg(:));
lon_bg=double(S.lon_bg(:));
v_bg=double(S.v_bg);
lon_obs=double(S.lon_obs(:));
observations=double(S.observations);
validate_input(time_bg,lon_bg,v_bg,lon_obs,observations,opt);

%% 2. Mooring-constrained weighted least-squares reconstruction
nlon=numel(lon_bg); ntime=numel(time_bg); nmoor=numel(lon_obs);
e=ones(nlon,1);
D2=spdiags([e -2*e e],[-1 0 1],nlon,nlon);
D2(1,:)=0; D2(end,:)=0;
H=zeros(nmoor,nlon);
for j=1:nmoor
    if lon_obs(j)<min(lon_bg) || lon_obs(j)>max(lon_bg)
        error('Mooring longitude %.2f is outside the GLORYS section.',lon_obs(j));
    end
    H(j,:)=interp1(lon_bg,eye(nlon),lon_obs(j),'linear',0);
end

lat0=58;
lon_start=-28;
xg=0:1:520;
lon_xg=lon_start+xg/(111.32*cosd(lat0));
v_fitted=nan(ntime,nlon);
v_on_1km=nan(ntime,numel(xg));
fit_residual_rms=nan(ntime,1);
A=[sqrt(opt.gamma)*H; ...
   sqrt(opt.lambda_background)*eye(nlon); ...
   sqrt(opt.lambda_smoothness)*D2];

for it=1:ntime
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
if ~any(all(isfinite(v_fitted),2))
    error('No complete assimilation profile was produced from the validation input.');
end

%% 3. Thirty-day quasi-steady background current
time=time_bg;
vg=v_on_1km;
dt_day=median(diff(time),'omitnan');
window_count=max(1,round(opt.background_window_days/dt_day));
if mod(window_count,2)==0, window_count=window_count+1; end
vg_smooth=movmean(vg,window_count,1,'omitnan');
time_day=(ceil(min(time)):floor(max(time)))';
vg_day=nan(numel(time_day),size(vg,2));
for j=1:numel(time_day)
    [~,nearest]=min(abs(time-time_day(j)));
    vg_day(j,:)=vg_smooth(nearest,:);
end

sample_note='';
if isfield(S,'sample_note'), sample_note=char(S.sample_note); end
source_description='';
if isfield(S,'source_description'), source_description=char(S.source_description); end
result=struct('latitude',lat0,'longitude',lon_bg, ...
    'background_velocity',v_bg,'background_time',time_bg, ...
    'mooring_longitude',lon_obs,'mooring_windowed_velocity',observations, ...
    'fitted_velocity',v_fitted,'fit_residual_rms',fit_residual_rms, ...
    'xg_km',xg,'lon_xg',lon_xg,'time',time,'vg',vg, ...
    'vg_smooth',vg_smooth,'time_day',time_day,'vg_day',vg_day, ...
    'observation_filter_window_days',14,'depth_range_m',[0 215], ...
    'options',opt,'input_file',input_file,'sample_note',sample_note, ...
    'source_description',source_description);

%% 4. Save and plot validation outputs
if opt.save_output
    save(fullfile(data_output_dir,'reanaly_fit_58N.mat'),'xg','lon_xg','vg','time', ...
        'vg_smooth','time_day','vg_day','window_count','opt','-v7');
end

if opt.make_figure
    iday=select_figure_day(time_bg,v_fitted,opt);
    fig=figure('Color','w','Units','inches','Position',[1 1 7.15 4.4]);
    plot(lon_bg,v_bg(iday,:),'k--','LineWidth',1.2); hold on;
    plot(lon_bg,v_fitted(iday,:),'b-','LineWidth',2);
    plot(lon_obs,observations(iday,:),'ro','MarkerSize',7,'LineWidth',1.3);
    yline(0,'Color',[0.65 0.65 0.65],'LineStyle',':');
    xlabel('Longitude (degree E)');
    ylabel('Meridional velocity (m s^{-1})');
    legend('GLORYS background','Mooring-constrained fit', ...
        'Mooring observations','Location','best');
    fit_date=datetime(time_bg(iday),'ConvertFrom','datenum','Format','yyyy-MM-dd');
    title(sprintf('58 N assimilation validation: %s',char(fit_date)));
    subtitle(sprintf('Fit residual RMS = %.4f m s^{-1}',fit_residual_rms(iday)));
    grid on; box off;
    result.figure_date=fit_date;
    result.figure_handle=fig;
    if opt.save_output
        exportgraphics(fig,fullfile(figure_output_dir,'reanalysis_fit_58N.png'), ...
            'Resolution',300);
    end
end
end

function validate_input(time_bg,lon_bg,v_bg,lon_obs,observations,opt)
if numel(time_bg)<3 || any(~isfinite(time_bg)) || any(diff(time_bg)<=0)
    error('time_bg must contain at least three finite, strictly increasing values.');
end
if numel(lon_bg)<3 || any(~isfinite(lon_bg)) || any(diff(lon_bg)<=0)
    error('lon_bg must contain at least three finite, strictly increasing values.');
end
if ~isequal(size(v_bg),[numel(time_bg),numel(lon_bg)])
    error('v_bg must have size [numel(time_bg), numel(lon_bg)].');
end
if ~isequal(size(observations),[numel(time_bg),numel(lon_obs)])
    error('observations must have size [numel(time_bg), numel(lon_obs)].');
end
if numel(time_bg)<opt.background_window_days
    warning('The sample is shorter than the requested background window.');
end
if any(~isfinite(lon_obs)), error('lon_obs must be finite.'); end
end

function iday=select_figure_day(time_bg,v_fitted,opt)
available=find(all(isfinite(v_fitted),2));
if isempty(available), error('No complete fitted profile is available for plotting.'); end
if ~isempty(opt.figure_date)
    if isdatetime(opt.figure_date)
        requested=datenum(opt.figure_date); %#ok<DATNM>
    elseif ischar(opt.figure_date) || isstring(opt.figure_date)
        requested=datenum(char(opt.figure_date)); %#ok<DATNM>
    elseif isnumeric(opt.figure_date) && isscalar(opt.figure_date)
        requested=double(opt.figure_date);
    else
        error('figure_date must be a datetime, date string, or MATLAB datenum.');
    end
    [~,q]=min(abs(time_bg(available)-requested));
    iday=available(q);
elseif ~isempty(opt.figure_index)
    requested=min(max(1,round(opt.figure_index)),numel(time_bg));
    [~,q]=min(abs(available-requested));
    iday=available(q);
else
    iday=available(ceil(numel(available)/2));
end
end

function opt=defaults(opt)
d=struct('gamma',400,'lambda_background',0.5,'lambda_smoothness',50, ...
    'background_window_days',30,'make_figure',true,'figure_date',[], ...
    'figure_index',[],'save_output',false);
names=fieldnames(d);
for j=1:numel(names)
    if ~isfield(opt,names{j}), opt.(names{j})=d.(names{j}); end
end
positive={'gamma','lambda_background','lambda_smoothness','background_window_days'};
for j=1:numel(positive)
    value=opt.(positive{j});
    if ~isscalar(value) || ~isfinite(value) || value<=0
        error('%s must be a positive finite scalar.',positive{j});
    end
end
if ~isempty(opt.figure_index) && (~isscalar(opt.figure_index) || ~isfinite(opt.figure_index))
    error('figure_index must be empty or one finite scalar.');
end
end
