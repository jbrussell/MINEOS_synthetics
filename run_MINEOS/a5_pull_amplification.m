% Pull values from .eig file for calculating amplification. Must first run
% a1 (mineos_nohang).
%
% !! THIS VERSION ONLY WORKS FOR FUNDAMENTAL MODE... 
% !! USE a5_pull_amplification_fromEig.m FOR OTHER MODE BRANCHES
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

clear

isfigure = 1;

%% get pamameters information 
parameter_FRECHET;
CARD = param.CARD;
CARDID =  param.CARDID;
periods = param.periods;

%% Set path to fortran executables
% setpath_mineos;


%% Change environment variables to deal with gfortran
setenv('GFORTRAN_STDIN_UNIT', '5') 
setenv('GFORTRAN_STDOUT_UNIT', '6') 
setenv('GFORTRAN_STDERR_UNIT', '0')


%% Run Spheroidal Branches First
if SONLY 
    
    % Write out run files -- be sure that paths are correct!
%     write_get_eigfxn_grvelo(periods);
    % Index surface layer (or seafloor)
    card=read_model_card([param.CARDPATH,'/',CARD]);
    ind = find(card.vsv~=0);
    ind_surf = ind(end);
    write_get_eigfxn_grvelo_int(periods,ind_surf)
    
    %mineos_nohang for s1-s5
    disp('Running get_eigfxn_grvelo');    
%     com = ['cat run_get_eigfxn_grvelo.s | get_eigfxn_grvelo'];
    com = ['cat run_get_eigfxn_grvelo.s | get_eigfxn_grvelo_int'];
    [status,log] = system(com);
% % %     fprintf(fid,log);
    if status ~= 0     
        error( 'something is wrong at get_eigfxn_grvelo')
    end
    
    delete('run_get_eigfxn_grvelo.s')
end

%% Run toroidal branches next
if TONLY
    error('Have not set this up for torroidal modes yet...');
end

% % Delete unnecessary files
% delete('*.LOG','qlog','log_table');
% delete([CARDTABLE,'log*']);
% delete([CARDTABLE,'*.asc']);
% if exist([CARDTABLE,param.CARDID,'.',TYPEID,'_0.eig_fix'],'file') ~= 0
%     delete([CARDTABLE,'*',TYPEID,'*.eig']); % Delete regular .eig files
% elseif num_loop == 0
%     delete([CARDTABLE,param.CARDID,'.',TYPEID,'_',num2str(num_loop),'.eig']);
% end

%% Change the environment variables back to the way they were
setenv('GFORTRAN_STDIN_UNIT', '-1') 
setenv('GFORTRAN_STDOUT_UNIT', '-1') 
setenv('GFORTRAN_STDERR_UNIT', '-1')
%% Calculate amplification

% AMP = load_eigfxn_grvelo_asc(CARDID,'S');
AMP = load_eigfxn_grvelo_int_asc(CARDID,'S');

% Eddy & Ekstrom (2014)
A_R = AMP.U_0 ./ sqrt(AMP.grv);
% Lin et al. (2012) and Bowden et al. (2017)
A_R_lin = (AMP.grv .* AMP.I_0).^(-0.5);

% % Calculate using eig file
% mbranch = 0; % mode number
% TYPE = 'S';
% eigmat = [param.eigpath,param.CARDID,'.',TYPE,num2str(mbranch),'.mat'];
% temp = load(eigmat);
% eig = temp.eig; clear temp;
% ind = find(eig(1).mu~=0);
% I_surf = ind(end);
% I_6371 = length(eig(1).z);
% for ip = 1:length(periods)
%     A_R_surf(ip) = eig(ip).u(I_surf) ./ sqrt(AMP.grv(ip));
%     
%     A_R_6371(ip) = eig(ip).u(I_6371) ./ sqrt(AMP.grv(ip));
% end

% % Calculate I_0 version
% for ip = 1:length(periods)
%     u = eig(ip).u;
%     v = eig(ip).v;
%     rho = eig(1).rho;
%     r = eig(1).r;
%     z = eig(1).z*1000;
%     
%     ind = find(eig(1).mu~=0);
%     I_surf = ind(end);
%     u = u ./ u(I_surf);
%     v = v ./ v(I_surf);
%     I_0 = trapz(r , rho.*((u).^2 + (v).^2).*r.^2);
% %     I_0 = trapz(r , rho.*((u).^2 + (v).^2));
% %     I_0 = abs(trapz(z , rho.*((u).^2 + (v).^2).*z.^2));
%     A_R_lin(ip) = (AMP.grv(ip) .* I_0).^(-0.5);
% end

if isfigure
    figure(99); clf;
    hold on; box on;
    plot(AMP.periods,A_R,'o-r','linewidth',2);
    % plot(periods,A_R_surf,'o-b','linewidth',2);
    % plot(periods,A_R_6371,'o--c','linewidth',2);
    plot(periods,A_R_lin,'o-m','linewidth',2);
    xlabel('Period (s)');
    ylabel('Amplification');
    set(gca,'fontsize',15,'linewidth',1.5);
    legend({'Eddy & Ekstrom (2014)','Lin et al. (2012)'},'location','eastoutside')
end






