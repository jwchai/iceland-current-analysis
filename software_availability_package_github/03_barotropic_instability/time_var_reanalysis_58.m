function result = time_var_reanalysis_58(velocity_file,topography_file,output_dir,options)
%TIME_VAR_REANALYSIS_58 Barotropic-instability diagnostics along 58 N.
%
% result = time_var_reanalysis_58(velocity_file,topography_file,output_dir,options)
%
% Required input variables
%   velocity_file: vg_smooth [time x x], time, xg
%   topography_file: topo_on_xg [1 x x]
%
% The default numerical configuration reproduces the reviewed annual test:
% a 1-km grid, a five-point topographic moving mean, f=1.23e-4 s^-1,
% wavelengths from 10 to 300 km in 10-km increments, and no post-solver
% velocity multiplier or modal-period subtraction. Following the manuscript
% Figure 4 criterion, a day is flagged when the maximum growth rate exceeds
% 0.1 d^-1 and the corresponding modal period is shorter than 20 days. The
% wavelength and e-folding propagation-distance ratio D/lambda are retained
% as diagnostics but do not determine the green event markers.
%
% Example
%   root = '/path/to/software_availability_package';
%   velocity_file = '/path/to/reanaly_fit_58N.mat';
%   topo_file = fullfile(root,'data','smoothed_topography_on_xg.mat');
%   out = fullfile(root,'outputs','barotropic_instability_58N');
%   result = time_var_reanalysis_58(velocity_file,topo_file,out,struct());

if nargin<1 || isempty(velocity_file)
    error('velocity_file is required and must be the output from module 01.');
end
if nargin<2 || isempty(topography_file)
    error('topography_file is required.');
end
if nargin<3, output_dir=''; end
if nargin<4 || isempty(options), options=struct(); end
opt = default_options(options);

validate_file(velocity_file,'velocity');
validate_file(topography_file,'topography');
if opt.save_output && isempty(output_dir)
    error('output_dir is required when save_output=true.');
end
if ~isempty(output_dir) && exist(output_dir,'dir')~=7, mkdir(output_dir); end

V = load(velocity_file);
T = load(topography_file);
require_fields(V,{'vg_smooth','time','xg'},'velocity file');
require_fields(T,{'topo_on_xg'},'topography file');

time = double(V.time(:));
velocity = double(V.vg_smooth);
xg_km = double(V.xg(:)');
topography = double(T.topo_on_xg(:)');
if size(velocity,1)~=numel(time) || size(velocity,2)~=numel(xg_km)
    error('vg_smooth must have size [numel(time), numel(xg)].');
end
if numel(topography)~=numel(xg_km)
    error('topo_on_xg and xg must have the same number of points.');
end
if any(diff(time)<=0), error('time must be finite and strictly increasing.'); end
if any(~isfinite(topography)), error('Topography input must be finite.'); end

start_dn = to_datenum_or_default(opt.start_date,-Inf);
end_dn = to_datenum_or_default(opt.end_date,Inf);
use = time>=start_dn & time<=end_dn;
time = time(use);
velocity = velocity(use,:);
if isempty(time), error('No velocity profiles fall within the requested dates.'); end
input_profile_valid=all(isfinite(velocity),2);
if ~any(input_profile_valid), error('No finite velocity profiles are available.'); end
if any(~input_profile_valid)
    warning('%d velocity profiles contain missing values and will return NaN diagnostics.', ...
        sum(~input_profile_valid));
end

dx_m = median(diff(xg_km))*1000;
if any(abs(diff(xg_km)*1000-dx_m)>1e-6)
    error('xg must be uniformly spaced.');
end
if abs(dx_m-1000)>1e-6
    error('The reviewed 58 N solver requires a 1-km xg grid.');
end

depth = -topography;
if any(depth<=0), error('-topo_on_xg must be positive water depth.'); end
depth = smoothdata(depth,'movmean',opt.topography_window_points);

nt = numel(time);
gr = nan(1,nt);
gr_l = nan(1,nt);
Wt = nan(1,nt);
last_completed = 0;
checkpoint_file = '';
if opt.save_output
    checkpoint_file = fullfile(output_dir,opt.output_filename);
    if opt.resume && exist(checkpoint_file,'file')==2
        P = load(checkpoint_file,'time','gr','gr_l','Wt','last_completed','options');
        if isfield(P,'options') && isequaln(P.time,time) && isequaln(P.options,opt)
            gr=P.gr; gr_l=P.gr_l; Wt=P.Wt; last_completed=P.last_completed;
            fprintf('Resuming after day %d of %d.\n',last_completed,nt);
        end
    end
end

use_parallel = opt.use_parallel && license('test','Distrib_Computing_Toolbox');
if opt.use_parallel && ~use_parallel
    warning('Parallel Computing Toolbox unavailable; using a serial loop.');
end
if use_parallel && isempty(gcp('nocreate')), parpool('local'); end

run_clock = tic;
for chunk_start=last_completed+1:opt.chunk_size:nt
    chunk_index=chunk_start:min(chunk_start+opt.chunk_size-1,nt);
    velocity_chunk=velocity(chunk_index,:);
    gr_chunk=nan(1,numel(chunk_index));
    wavelength_chunk=nan(1,numel(chunk_index));
    period_chunk=nan(1,numel(chunk_index));
    if use_parallel
        parfor j=1:numel(chunk_index)
            [gr_chunk(j),wavelength_chunk(j),period_chunk(j)] = ...
                diagnose_profile(velocity_chunk(j,:),depth,dx_m,opt);
        end
    else
        for j=1:numel(chunk_index)
            [gr_chunk(j),wavelength_chunk(j),period_chunk(j)] = ...
                diagnose_profile(velocity_chunk(j,:),depth,dx_m,opt);
        end
    end
    gr(chunk_index)=gr_chunk;
    gr_l(chunk_index)=wavelength_chunk;
    Wt(chunk_index)=period_chunk;
    last_completed=chunk_index(end);
    if opt.save_output
        options=opt;
        save(checkpoint_file,'time','gr','gr_l','Wt','last_completed','options', ...
            'velocity_file','topography_file');
    end
    elapsed=toc(run_clock);
    eta=elapsed/max(1,last_completed)*(nt-last_completed);
    fprintf('Completed %d/%d days (%.1f%%); elapsed %.1f min; ETA %.1f min.\n', ...
        last_completed,nt,100*last_completed/nt,elapsed/60,eta/60);
end

D_over_lambda=1./(gr.*Wt);
D_over_lambda(~isfinite(D_over_lambda) | gr<=0 | Wt<=0)=NaN;
unstable=gr>opt.growth_threshold_day & ...
    Wt>0 & Wt<opt.modal_period_threshold_day;
locally_amplifying=D_over_lambda<opt.distance_ratio_threshold;
time_cond=time(unstable);

result=struct('time',time,'growth_rate_day',gr,'wavelength_km',gr_l, ...
    'earth_fixed_period_day',Wt,'efolding_distance_over_wavelength',D_over_lambda, ...
    'unstable',unstable,'locally_amplifying',locally_amplifying, ...
    'time_cond',time_cond,'options',opt, ...
    'velocity_file',velocity_file,'topography_file',topography_file, ...
    'input_profile_valid',input_profile_valid,'last_completed',last_completed);

if opt.save_output
    save(checkpoint_file,'time','gr','gr_l','Wt','D_over_lambda','unstable', ...
        'locally_amplifying','time_cond','input_profile_valid','last_completed','options', ...
        'velocity_file','topography_file');
end
if opt.make_figure
    plot_diagnostics(result,output_dir,opt.save_output);
end
end

function [max_growth,wavelength,period] = diagnose_profile(vm,depth,dx,opt)
if any(~isfinite(vm))
    max_growth=NaN; wavelength=NaN; period=NaN; return;
end
dvm=gradient(vm)/dx;
dh=gradient(depth)/dx;
[modes,wavelengths]=solve_modes(vm,depth,dx,opt.coriolis_s,opt.gravity_m_s2, ...
    opt.wavelength_min_km,opt.wavelength_max_km,opt.wavelength_step_km,dh,dvm);
growth=squeeze(modes(:,5,:));
[max_growth,linear_index]=max(growth(:),[],'omitnan');
if isempty(linear_index) || ~isfinite(max_growth)
    max_growth=NaN; wavelength=NaN; period=NaN; return;
end
[mode_index,wavelength_index]=ind2sub(size(growth),linear_index);
wavelength=wavelengths(wavelength_index);
period=modes(mode_index,1,wavelength_index);
end

function [output,wavelengths] = solve_modes(vm,h,dx,f,g,L1,L2,DL,dh,dvm)
nmax=numel(vm);
ndim=3*nmax-3;
wavelengths=L2:-DL:L1;
output=nan(ndim,5,numel(wavelengths));

for il=1:numel(wavelengths)
    LL=wavelengths(il);
    k=2*pi/(LL*1000);
    A=zeros(ndim);

    % Free-surface equation.
    for j=1:nmax-1
        A(j,j)=k*vm(j);
        if j==1
            A(nmax,j)=h(1)/dx-0.5*dh(1);
        elseif j<nmax-1
            A(j+nmax-1,j)=h(j)/dx-0.5*dh(j);
            A(j+nmax,j)=-h(j)/dx-0.5*dh(j);
        end
        A(2*nmax-2+j,j)=h(j)*k/2;
        if 2*nmax-1+j<=ndim, A(2*nmax-1+j,j)=h(j)*k/2; end
    end
    A(nmax+1,1)=-h(1)/dx-0.5*dh(1);
    A(2*nmax-1,nmax-1)=0;
    A(2*nmax-2,nmax-1)=h(nmax-1)/dx-0.5*dh(nmax-1);

    % Along-section momentum equation.
    for j=nmax:2*nmax-2
        i=j-nmax;
        if i>=1
            A(i+1,j)=g/dx;
            A(i,j)=-g/dx;
        end
        A(j,j)=k*vm(j-nmax+1);
        target=j+nmax-1;
        if target<=ndim, A(target,j)=-f; end
    end

    % Cross-section momentum equation.
    for j=2*nmax-1:ndim
        i=j-2*nmax+1;
        if i>=1
            A(i,j-1)=g*k/2;
            A(i,j)=g*k/2;
        end
        A(nmax-1,ndim)=g*k;
        target=j-nmax+1;
        if target>=nmax+1 && target<=2*nmax-2
            A(target,j)=-f-dvm(target-nmax);
        end
        A(j,j)=k*vm(j-2*nmax+2);
    end
    A(nmax,2*nmax-1)=-f-dvm(1);

    eigenvalue=eig(A);
    angular_frequency=abs(real(eigenvalue));
    growth_frequency=abs(imag(eigenvalue));
    period_hour=2*pi./angular_frequency/3600;
    growth_period_hour=2*pi./growth_frequency/3600;
    inertial_period_hour=2*pi/f/3600;
    selected=period_hour>=inertial_period_hour & period_hour<=180*24;
    count=sum(selected);
    output(1:count,1,il)=period_hour(selected)/24;
    output(1:count,2,il)=growth_period_hour(selected)/24;
    output(1:count,3,il)=angular_frequency(selected)/k;
    output(1:count,4,il)=LL./output(1:count,3,il)/3.6;
    output(1:count,5,il)=imag(eigenvalue(selected))*86400;
end
end

function plot_diagnostics(result,output_dir,save_output)
date=datetime(result.time,'ConvertFrom','datenum');
fig=figure('Color','w','Units','inches','Position',[1 1 7.15 6.2]);
tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
labels={'Growth rate (d^{-1})','Wavelength (km)','Earth-fixed period (days)'};
titles={'a  Maximum modal growth rate','b  Wavelength of the most unstable mode', ...
    'c  Period of the most unstable mode'};
values={result.growth_rate_day,result.wavelength_km,result.earth_fixed_period_day};
colors={[0.835 0.369 0],[0 0.447 0.698],[0.25 0.25 0.25]};
for panel=1:3
    ax=nexttile; hold(ax,'on');
    plot(ax,date,values{panel},'Color',colors{panel},'LineWidth',1.2);
    scatter(ax,date(result.unstable),values{panel}(result.unstable),8, ...
        [0 0.62 0.451],'filled','MarkerFaceAlpha',0.7);
    ylabel(ax,labels{panel}); title(ax,titles{panel},'FontWeight','normal');
    grid(ax,'on'); box(ax,'off');
end
xlabel('Time');
if save_output
    exportgraphics(fig,fullfile(output_dir,'barotropic_instability_58N.png'),'Resolution',300);
    exportgraphics(fig,fullfile(output_dir,'barotropic_instability_58N.pdf'),'ContentType','vector');
end
end

function opt=default_options(opt)
d=struct('coriolis_s',1.23e-4,'gravity_m_s2',9.81, ...
    'wavelength_min_km',10,'wavelength_max_km',300,'wavelength_step_km',10, ...
    'topography_window_points',5,'growth_threshold_day',0.1, ...
    'modal_period_threshold_day',20,'distance_ratio_threshold',0.25, ...
    'start_date',[],'end_date',[],'chunk_size',30,'use_parallel',true, ...
    'resume',true,'make_figure',true,'save_output',true, ...
    'output_filename','barotropic_instability_58N.mat');
names=fieldnames(d);
for j=1:numel(names)
    if ~isfield(opt,names{j}) || isempty(opt.(names{j})), opt.(names{j})=d.(names{j}); end
end
if opt.wavelength_max_km<=opt.wavelength_min_km || opt.wavelength_step_km<=0
    error('Wavelength bounds and step are invalid.');
end
if mod(opt.wavelength_max_km-opt.wavelength_min_km,opt.wavelength_step_km)~=0
    error('Wavelength range must be divisible by wavelength_step_km.');
end
if ~isscalar(opt.modal_period_threshold_day) || opt.modal_period_threshold_day<=0
    error('modal_period_threshold_day must be a positive scalar.');
end
end

function value=to_datenum_or_default(input,default_value)
if isempty(input), value=default_value;
elseif isdatetime(input), value=datenum(input); %#ok<DATNM>
else, value=double(input);
end
end

function validate_file(file,label)
if exist(file,'file')~=2, error('%s file not found: %s',label,file); end
end

function require_fields(S,names,label)
for j=1:numel(names)
    if ~isfield(S,names{j}), error('%s is missing variable %s.',label,names{j}); end
end
end
