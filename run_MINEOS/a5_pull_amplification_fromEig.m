% Pull values from .eig file for calculating amplification. Must first run
% a1 (mineos_nohang) and a4 (to generate eigenfunctions). The advantage of
% this version of "pull_amplification" is that it can be used for mode
% branches other than the fundamental mode.
%
% This includes two different parameterizations of amplification
% 1) Eddy & Ekstrom (2014)
%    A_R = U_0 / sqrt(grv)
%
% 2) Lin et al. (2012)
%    A_R = (grv * I_0)^-0.5
%          where
%        I_0 = integral( rho * (U^2 + V^2) * r ) dr
%
% Their values appear quite different, but they are nearly identical when 
% compared to a reference station A_R0. 
%
% In other words, the quantity
%    A_R / A_R0    (where A_R0 is some arbitrary reference profile)
% is nearly the same for both parameterizations
%
% jbrussell - 5/4/2022
%

clear

mbranch = 0; % mode branch of interest

isfigure = 1;

%% get pamameters information 
parameter_FRECHET;
CARD = param.CARD;
CARDID =  param.CARDID;

%% Run Spheroidal Branches First
if SONLY 
    TYPEID = param.STYPEID;
    TYPE = 'S';
    eigmat = [param.eigpath,param.CARDID,'.',TYPE,num2str(mbranch),'.mat'];
    temp = load(eigmat);
    eig = temp.eig; clear temp;
end

%% Run toroidal branches next
if TONLY
    error('Have not set this up for torroidal modes yet...');
    TYPEID = param.TTYPEID;
    TYPE = 'T';
    eigmat = [param.eigpath,param.CARDID,'.',TYPE,num2str(mbranch),'.mat'];
    temp = load(eigmat);
    eig = temp.eig; clear temp;
end

%% Calculate amplification

periods = [eig.per_want];
qfile = [param.TABLEPATH,param.CARDID,'/tables/',param.CARDID,'.',TYPEID,'.q'];
[phv,grv,phvq] = readMINEOS_qfile_per(qfile,periods,mbranch);

% Eddy & Ekstrom (2014)
ind = find(eig(1).mu~=0);
I_surf = ind(end);
I_6371 = length(eig(1).z);
for ip = 1:length(periods)
    A_R_surf(ip) = eig(ip).u(I_surf) ./ sqrt(grv(ip));
    
%     A_R_6371(ip) = eig(ip).u(I_6371) ./ sqrt(grv(ip));
end

% Calculate I_0 version
for ip = 1:length(periods)
    u = eig(ip).u;
    v = eig(ip).v;
    rho = eig(1).rho;
    r = eig(1).r;
    z = eig(1).z*1000;
    
    ind = find(eig(1).mu~=0);
    I_surf = ind(end);
    u = u ./ u(I_surf);
    v = v ./ v(I_surf);
    I_0 = trapz(r , rho.*((u).^2 + (v).^2).*r.^2);
%     I_0 = trapz(r , rho.*((u).^2 + (v).^2));
%     I_0 = abs(trapz(z , rho.*((u).^2 + (v).^2).*z.^2));
    A_R_lin(ip) = (grv(ip) .* I_0).^(-0.5);
end

if isfigure
    figure(99); clf;
    hold on; box on;
%     plot(periods,A_R,'o-r','linewidth',2);
    plot(periods,A_R_surf,'o-b','linewidth',2);
%     plot(periods,A_R_6371,'o--c','linewidth',2);
    plot(periods,A_R_lin,'o-m','linewidth',2);
    xlabel('Period (s)');
    ylabel('Amplification');
    set(gca,'fontsize',15,'linewidth',1.5);
    legend({'Eddy & Ekstrom (2014)','Lin et al. (2012)'},'location','eastoutside')
end






