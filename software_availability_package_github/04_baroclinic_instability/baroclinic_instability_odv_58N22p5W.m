function result=baroclinic_instability_odv_58N22p5W(input_file,options)
%BAROCLINIC_INSTABILITY_ODV_58N22P5W Package icelandcurrent.m lines 103--224.

if nargin<1 || isempty(input_file)
    module_dir=fileparts(mfilename('fullpath'));
    input_file=fullfile(module_dir,'data','odv_58N_22p5W_validation.mat');
end
if nargin<2 || isempty(options), options=struct(); end
opt=default_options(options);

module_dir=fileparts(mfilename('fullpath'));
data_dir=fullfile(module_dir,'data');
figure_dir=fullfile(module_dir,'figure');
validate_dependencies();
if exist(input_file,'file')~=2, error('Missing input: %s',input_file); end
D=load(input_file);

%% Potential density
t=double(D.temperature_degC(:));
s=double(D.salinity_psu(:));
p=double(D.pressure_source_dbar(:));
pden=sw_pden(s,t,p,0.4e2);

%% Thermal-wind velocity
depth=double(D.depth_m(:));
sigma0=double(D.sigma0_stencil);
depth(depth==0)=NaN;
sigma0(sigma0==0)=NaN;

dist_x=1e3*sw_dist([40,40],[300,300.5],'km');
dist_y=1e3*sw_dist([40,40.5],[300,300],'km');
f=sw_f(58);
depth=-depth;
for iz=1:size(sigma0,3)
    [sigma0y_layer,sigma0x_layer]=gradient(sigma0(:,:,iz));
    sigma0x(iz)=sigma0x_layer(2,2)/dist_x; %#ok<AGROW>
    sigma0y(iz)=sigma0y_layer(2,2)/dist_y; %#ok<AGROW>
end
sigma0x=sigma0x(1:38);
sigma0y=sigma0y(1:38);
depth=depth(1:38);
depth(1)=0;

for iz=23:38
    vg(iz)=-1*9.8*trapz(depth(22:iz),sigma0x(22:iz))/f*1e-3; %#ok<AGROW>
    ug(iz)= 1*9.8*trapz(depth(22:iz),sigma0y(22:iz))/f*1e-3; %#ok<AGROW>
end
depth1=flip(depth(1:22));
sigma0x1=flip(sigma0x(1:22));
sigma0y1=flip(sigma0y(1:22));
for iz=2:22
    vg1(iz)=-1*9.8*trapz(depth1(1:iz),sigma0x1(1:iz))/f*1e-3; %#ok<AGROW>
    ug1(iz)= 1*9.8*trapz(depth1(1:iz),sigma0y1(1:iz))/f*1e-3; %#ok<AGROW>
end
vg(1:22)=flip(vg1);
ug(1:22)=flip(ug1);

%% Instability analysis
f0=sw_f(58);
rho=pden;
rho0=1024;
g=9.81;
z=depth';
U=ug;
V=vg;
F=f0^2*rho0/g;
betax=0;
betay=(sw_f(40.5)-sw_f(40))/sw_dist([39.5,40],[1,1],'km')*1e-3;
kvec=linspace(-2e-3,2e-3,100);
lvec=linspace(-2e-3,2e-3,100);
[wi_max,wr_max,psiv,amp,phase]=qggrz_source( ...
    z(1:31),rho(1:31),U(1:31),V(1:31),F,betax,betay,kvec,lvec,1);

n=0;
locx=[];
locy=[];
for ii=2:size(wi_max,1)-1
    for jj=2:size(wi_max,2)-1
        if wi_max(ii,jj)>wi_max(ii-1,jj) && ...
                wi_max(ii,jj)>wi_max(ii,jj-1) && ...
                wi_max(ii,jj)>wi_max(ii+1,jj) && ...
                wi_max(ii,jj)>wi_max(ii,jj+1)
            n=n+1;
            locx(n)=ii; %#ok<AGROW>
            locy(n)=jj; %#ok<AGROW>
        end
    end
end

result=struct('longitude',double(D.longitude),'latitude',double(D.latitude), ...
    'source_grid_indices_matlab',double(D.source_grid_indices_matlab), ...
    'pden',pden,'depth',depth,'sigma0x',sigma0x,'sigma0y',sigma0y, ...
    'ug',ug,'vg',vg,'f0',f0,'F',F,'betax',betax,'betay',betay, ...
    'kvec',kvec,'lvec',lvec,'wi_max',wi_max,'wr_max',wr_max, ...
    'psiv',psiv,'amp',amp,'phase',phase,'locx',locx,'locy',locy);

if opt.save_output
    if exist(data_dir,'dir')~=7, mkdir(data_dir); end
    save(fullfile(data_dir,opt.output_filename),'result','-v7.3');
end
if opt.make_figure
    if exist(figure_dir,'dir')~=7, mkdir(figure_dir); end
    plot_growth(result,figure_dir,opt);
end
end

function [wi_max,wr_max,psiv,amp,phase]= ...
    qggrz_source(z,rho,U,V,F,betax,betay,kvec,lvec,dim)
% The ten-input qggrz formulation called by icelandcurrent.m.
nkx=length(kvec);
nky=length(lvec);
z=z(:);
rho=rho(:);
U=U(:);
V=V(:);
nz=max([length(U),length(V),length(z),length(rho)]);
G=pv_stretch_opz_source(z,rho,F,dim);
Q_x=betax+G*V;
Q_y=betay-G*U;
wi_max=zeros(nkx,nky);
wr_max=zeros(nkx,nky);
psiv=complex(zeros(nkx,nky,nz));
amp=zeros(nkx,nky,nz);
phase=zeros(nkx,nky,nz);
for kc=1:nkx
    k=kvec(kc);
    for lc=1:nky
        l=lvec(lc);
        K2=k^2+l^2;
        KdotU=k*U+l*V;
        KxdelQ=k*Q_y-l*Q_x;
        B=G-K2*eye(nz);
        A=diag(KxdelQ)+diag(KdotU)*B;
        [evec,D]=eig(A,B);
        w=diag(D);
        if any(imag(w)>0)
            [wi_max(kc,lc),ind]=max(imag(w));
            wr_max(kc,lc)=real(w(ind));
        else
            [wr_max(kc,lc),ind]=max(real(w));
            wi_max(kc,lc)=imag(w(ind));
        end
        psiv(kc,lc,:)=evec(:,ind);
        amp(kc,lc,:)=sqrt(real(evec(:,ind)).^2+imag(evec(:,ind)).^2);
        phase(kc,lc,:)=atan2(imag(evec(:,ind)),real(evec(:,ind)));
    end
end
end

function G=pv_stretch_opz_source(z,rho,F,dim)
nz=length(z);
rho=rho(:);
z=z(:);
Dz=z(1:end-1)-z(2:end);
dz=[Dz(1);(Dz(1:end-1)+Dz(2:end))/2;Dz(end)];
drho=rho(2:end)-rho(1:end-1);
if dim==0
    drho0=(rho(end)-rho(1))/(nz-1);
    drho=drho/drho0;
    dz=dz/sum(dz);
end
sub=F./(dz(2:end).*drho);
sup=F./(dz(1:end-1).*drho);
mid=zeros(nz,1);
mid(1)=-F/(dz(1)*drho(1));
sup(1)=F/(dz(1)*drho(1));
mid(end)=-F/(dz(end)*drho(end));
sub(end)=F/(dz(end)*drho(end));
if nz>2
    mid(2:end-1)=-F./dz(2:end-1).* ...
        (1./drho(1:end-1)+1./drho(2:end));
end
G=diag(sub,-1)+diag(mid)+diag(sup,1);
end

function plot_growth(result,figure_dir,opt)
fig=figure('Color','w','Visible',opt.figure_visible);
[~,h1]=contourf(result.kvec*1e3,result.lvec*1e3, ...
    result.wi_max*3600*24,200);
h1.LineStyle='none';
colorbar;
colormap(flip(othercolor('Spectral11')));
caxis([0 0.04]);
set(gca,'LineWidth',2,'FontSize',16);
xlabel('k [1/km]');
ylabel('l [1/km]');
hold on;
scatter(result.kvec(result.locy([2,3,4]))*1e3, ...
    result.lvec(result.locx([2,3,4]))*1e3,80,'k','filled','^');
scatter(result.kvec(result.locy([9,10,11]))*1e3, ...
    result.lvec(result.locx([9,10,11]))*1e3,80,'k','filled','^');
xticks([-2 -1 0 1 2]);
yticks([-2 -1 0 1 2]);
exportgraphics(fig,fullfile(figure_dir,[opt.figure_basename '.png']), ...
    'Resolution',300);
if opt.save_pdf
    exportgraphics(fig,fullfile(figure_dir,[opt.figure_basename '.pdf']), ...
        'ContentType','vector');
end
if strcmpi(opt.figure_visible,'off'), close(fig); end
end

function opt=default_options(opt)
defaults=struct('save_output',true,'make_figure',true, ...
    'output_filename','baroclinic_instability_odv_validation.mat', ...
    'figure_basename','baroclinic_instability_odv_validation', ...
    'figure_visible','off','save_pdf',true);
names=fieldnames(defaults);
for i=1:numel(names)
    if ~isfield(opt,names{i}) || isempty(opt.(names{i}))
        opt.(names{i})=defaults.(names{i});
    end
end
end

function validate_dependencies()
required={'sw_pden','sw_f','sw_dist','othercolor'};
missing=required(cellfun(@(name) exist(name,'file')~=2,required));
if ~isempty(missing)
    error('Missing source-program functions: %s',strjoin(missing,', '));
end
end
