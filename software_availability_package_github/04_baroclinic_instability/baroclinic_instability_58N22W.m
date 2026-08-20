function result = baroclinic_instability_58N22W(time_index,options)
%BAROCLINIC_INSTABILITY_58N22W Original 58 N, 22 W QG workflow.
%
% result = baroclinic_instability_58N22W(time_index,options)
%
% This public function preserves the analysis sequence of the author's
% build_instability_dataset_58N22W.m and qggrz.m programs. Input and output
% paths are resolved relative to this file:
%
%   04_baroclinic_instability/
%     baroclinic_instability_58N22W.m
%     data/uvts.nc
%     data/glo_topo.nc
%     outputs/
%
% The default time_index=16 reproduces the case selected in the original
% script. Set time_index to another scalar index to diagnose another monthly
% mean. Required dependency: MATLAB Seawater Toolbox (sw_f, sw_dist,
% sw_pres, sw_bfrq, and sw_pden).

if nargin<1 || isempty(time_index), time_index=16; end
if nargin<2 || isempty(options), options=struct(); end
opt=default_options(options);
if ~isscalar(time_index) || ~isfinite(time_index) || time_index~=round(time_index)
    error('time_index must be one finite integer.');
end

script_dir=fileparts(mfilename('fullpath'));
data_dir=fullfile(script_dir,'data');
output_dir=fullfile(script_dir,'outputs');
uvfile=fullfile(data_dir,'uvts.nc');
topofile=fullfile(data_dir,'glo_topo.nc');
validate_inputs(uvfile,topofile);
if (opt.save_output || opt.make_figure) && exist(output_dir,'dir')~=7
    mkdir(output_dir);
end

%% Coordinates and selected monthly-mean column
lon=double(ncread(uvfile,'longitude')); lon=lon(:);
lat=double(ncread(uvfile,'latitude')); lat=lat(:);
z=double(ncread(uvfile,'depth')); z=z(:);
t=double(ncread(uvfile,'time')); t=t(:);
nx=numel(lon); ny=numel(lat); nz=numel(z); nt=numel(t);
if time_index<1 || time_index>nt
    error('time_index=%d is outside the available range 1:%d.',time_index,nt);
end
[~,ii]=min(abs(lon-opt.longitude));
[~,jj]=min(abs(lat-opt.latitude));
if ii<=1 || ii>=nx || jj<=1 || jj>=ny
    error('Target point is on the data-grid edge.');
end
lon_pt=lon(ii); lat_pt=lat(jj);

thetao=read_column(uvfile,'thetao',ii,jj,nz,time_index);
so=read_column(uvfile,'so',ii,jj,nz,time_index);
uo=read_column(uvfile,'uo',ii,jj,nz,time_index);
vo=read_column(uvfile,'vo',ii,jj,nz,time_index);

%% Local rotation, beta, bathymetry, and bottom slope
f0=sw_f(lat_pt);
earth_rotation=7.292115e-5;
earth_radius=6371000;
beta=2*earth_rotation*cosd(lat_pt)/earth_radius;

deptho=double(ncread(topofile,'deptho'));
if isequal(size(deptho),[ny nx]), deptho=deptho'; end
if ~isequal(size(deptho),[nx ny])
    error('deptho dimensions do not match the longitude-latitude grid.');
end
dx_w=1000*first_value(sw_dist([lat_pt lat_pt],[lon(ii-1) lon_pt],'km'));
dx_e=1000*first_value(sw_dist([lat_pt lat_pt],[lon_pt lon(ii+1)],'km'));
dy_s=1000*first_value(sw_dist([lat(jj-1) lat_pt],[lon_pt lon_pt],'km'));
dy_n=1000*first_value(sw_dist([lat_pt lat(jj+1)],[lon_pt lon_pt],'km'));
bottom_height=-deptho;
bottom_hx=(bottom_height(ii+1,jj)-bottom_height(ii-1,jj))/(dx_e+dx_w);
bottom_hy=(bottom_height(ii,jj+1)-bottom_height(ii,jj-1))/(dy_n+dy_s);

%% Density, stratification, and valid full-depth profiles
pressure=sw_pres(z,lat_pt);
density=double(sw_pden(so,thetao,pressure,0)); density=density(:);
[N2_mid,~]=sw_bfrq(so,thetao,pressure,lat_pt); N2_mid=double(N2_mid(:));
valid=isfinite(z)&isfinite(density)&isfinite(uo)&isfinite(vo);
z_used=z(valid); density_used=density(valid);
u_used=uo(valid); v_used=vo(valid);
if numel(z_used)<5, error('Fewer than five valid vertical levels remain.'); end
if any(diff(z_used)<=0), error('Depth must increase downward.'); end
if any(diff(density_used)<=0)
    error('Selected potential-density profile is not strictly stable.');
end
% The gridded column terminates below the local seafloor, so its deepest
% nominal N2 entries are fill values. Average the deepest finite 5% rather
% than the final array slots (which made the author's scalar Fb become NaN).
N2_finite=N2_mid(isfinite(N2_mid)&N2_mid>0);
n_bottom=max(2,round(opt.bottom_fraction*numel(N2_finite)));
if numel(N2_finite)<n_bottom, error('Insufficient positive finite N2 values.'); end
N2_bottom=mean(N2_finite(end-n_bottom+1:end));
bottom=struct('hx',bottom_hx,'hy',bottom_hy,'N2b',N2_bottom,'f0',f0);

%% Author's continuously stratified QG eigenvalue calculation
F=f0^2*opt.reference_density/opt.gravity;
[growth,frequency,eigenvector,amplitude,phase,residual]=qggrz_original( ...
    -z_used,density_used,u_used,v_used,F,0,beta,opt.kvec,opt.lvec,bottom);

[maximum_growth_s,index]=max(growth(:));
[ik,il]=ind2sub(size(growth),index);
fastest_k=opt.kvec(ik); fastest_l=opt.lvec(il);
total_wavenumber=hypot(fastest_k,fastest_l);
if maximum_growth_s>0 && total_wavenumber>0
    inverse_wavenumber_scale_km=1/total_wavenumber/1000;
    wavelength_km=2*pi/total_wavenumber/1000;
else
    inverse_wavenumber_scale_km=NaN; wavelength_km=NaN;
end

result=struct('time_index',time_index,'time_value',t(time_index), ...
    'longitude',lon_pt,'latitude',lat_pt,'z',z_used,'density',density_used, ...
    'u',u_used,'v',v_used,'f0',f0,'beta',beta,'bottom',bottom, ...
    'kvec',opt.kvec,'lvec',opt.lvec,'growth_s',growth, ...
    'growth_per_day',growth*86400,'frequency_s',frequency, ...
    'eigenvector',eigenvector,'amplitude',amplitude,'phase',phase, ...
    'relative_residual',residual,'maximum_growth_per_day',maximum_growth_s*86400, ...
    'fastest_k',fastest_k,'fastest_l',fastest_l, ...
    'inverse_wavenumber_scale_km',inverse_wavenumber_scale_km, ...
    'wavelength_km',wavelength_km,'options',opt, ...
    'uv_file_relative',fullfile('data','uvts.nc'), ...
    'topography_file_relative',fullfile('data','glo_topo.nc'));

if opt.save_output
    save(fullfile(output_dir,sprintf('baroclinic_instability_t%02d.mat',time_index)), ...
        'result','-v7.3');
end
if opt.make_figure
    plot_growth(result,output_dir,opt.save_output);
end
end

function profile=read_column(file,name,ii,jj,nz,it)
profile=squeeze(double(ncread(file,name,[ii jj 1 it],[1 1 nz 1])));
profile=profile(:);
end

function value=first_value(value)
value=value(1);
end

function [growth,frequency,eigenvector,amplitude,phase,residual_max]= ...
    qggrz_original(z,rho,U,V,F,beta_x,beta_y,kvec,lvec,bottom)
% Original qggrz formulation with the accidental duplicate boundary update
% removed. The optional topographic bottom shear is passed explicitly.
z=z(:); rho=rho(:); U=U(:); V=V(:); nz=numel(z);
G=pv_stretch_opz_original(z,rho,F);
Q_x=beta_x+G*V;
Q_y=beta_y-G*U;

dzv=diff(z);
delta=abs([dzv(1);0.5*(dzv(1:end-1)+dzv(2:end));dzv(end)]);
drho=diff(rho);
alpha_top=F/(delta(1)*drho(1));
alpha_bottom=F/(delta(end)*drho(end));
dU_top=(U(2)-U(1))/delta(1);
dV_top=(V(2)-V(1))/delta(1);
dU_bottom=(U(end)-U(end-1))/delta(end)+(bottom.N2b/bottom.f0)*bottom.hy;
dV_bottom=(V(end)-V(end-1))/delta(end)-(bottom.N2b/bottom.f0)*bottom.hx;

nk=numel(kvec); nl=numel(lvec);
growth=zeros(nk,nl); frequency=zeros(nk,nl);
eigenvector=complex(nan(nk,nl,nz)); residual_max=nan(nk,nl);
for ik=1:nk
    k=kvec(ik);
    for il=1:nl
        l=lvec(il); K2=k^2+l^2;
        if K2==0, continue; end
        B=G-K2*eye(nz);
        A=diag(k*Q_y-l*Q_x)+diag(k*U+l*V)*B;
        % Bretherton/Smith boundary delta-PV terms: exactly once.
        A(1,1)=A(1,1)+alpha_top*(k*dV_top-l*dU_top);
        A(end,end)=A(end,end)+alpha_bottom*(k*dV_bottom-l*dU_bottom);
        [modes,D]=eig(A,B); omega=diag(D);
        mode_norm=vecnorm(modes,2,1)';
        denominator=(norm(A,1)+abs(omega)*norm(B,1)).*mode_norm;
        residual=vecnorm(A*modes-(B*modes).*omega.',2,1)'./max(denominator,eps);
        acceptable=isfinite(omega)&residual<1e-8;
        if any(acceptable)
            ids=find(acceptable); [~,q]=max(imag(omega(ids))); selected=ids(q);
            growth(ik,il)=max(0,imag(omega(selected)));
            frequency(ik,il)=real(omega(selected));
            eigenvector(ik,il,:)=modes(:,selected);
            residual_max(ik,il)=residual(selected);
        end
    end
end
amplitude=abs(eigenvector); phase=angle(eigenvector);
end

function G=pv_stretch_opz_original(z,rho,F)
% Author's dimensional stretching operator used by qggrz.
nz=numel(z); dz=get_dz_original(z); drho=diff(rho);
if any(dz<=0) || any(drho<=0)
    error('The original stretching operator requires positive layer thickness and density increments.');
end
sub=F./(dz(2:end).*drho);
sup=F./(dz(1:end-1).*drho);
mid=zeros(nz,1);
mid(1)=-F/(dz(1)*drho(1)); sup(1)=F/(dz(1)*drho(1));
mid(end)=-F/(dz(end)*drho(end)); sub(end)=F/(dz(end)*drho(end));
if nz>2
    mid(2:end-1)=-F./dz(2:end-1).*(1./drho(1:end-1)+1./drho(2:end));
end
G=diag(sub,-1)+diag(mid)+diag(sup,1);
end

function dz=get_dz_original(z)
Dz=z(1:end-1)-z(2:end);
dz=[Dz(1);0.5*(Dz(1:end-1)+Dz(2:end));Dz(end)];
end

function plot_growth(result,output_dir,save_output)
figure('Color','w');
contourf(result.kvec*1e3,result.lvec*1e3,result.growth_per_day',200, ...
    'LineStyle','none');
axis square; colorbar; colormap(turbo);
xlabel('k (km^{-1})'); ylabel('l (km^{-1})');
title(sprintf('Baroclinic growth rate, time index %d (d^{-1})',result.time_index));
if save_output
    exportgraphics(gcf,fullfile(output_dir,sprintf('baroclinic_growth_t%02d.png', ...
        result.time_index)),'Resolution',600);
end
end

function opt=default_options(opt)
d=struct('longitude',-22,'latitude',58,'gravity',9.81,'reference_density',1024, ...
    'bottom_fraction',0.05,'kvec',linspace(-2e-5,2e-5,100), ...
    'lvec',linspace(-2e-5,2e-5,100),'make_figure',true,'save_output',true);
names=fieldnames(d);
for j=1:numel(names)
    if ~isfield(opt,names{j}) || isempty(opt.(names{j})), opt.(names{j})=d.(names{j}); end
end
opt.kvec=double(opt.kvec(:)'); opt.lvec=double(opt.lvec(:)');
end

function validate_inputs(uvfile,topofile)
if exist(uvfile,'file')~=2, error('Missing relative input: %s',uvfile); end
if exist(topofile,'file')~=2, error('Missing relative input: %s',topofile); end
required={'sw_f','sw_dist','sw_pres','sw_bfrq','sw_pden'};
for j=1:numel(required)
    if exist(required{j},'file')~=2
        error('%s from the MATLAB Seawater Toolbox is required.',required{j});
    end
end
end
