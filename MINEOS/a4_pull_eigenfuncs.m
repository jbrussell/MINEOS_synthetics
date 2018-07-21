%% Plot eigenfunctions
% 9/14/15: JOSH RUSSELL
% 
% This program loads and plots the eigenfunctions
%
% MUST FIRST RUN 'run_mineos.m' to generate eigenfunction file and
% 'mk_kernels.m' to generate branch file
%
%
% calls on programs: 
%           parameter_FRECHET
%           eigen_ascii
%           
% JBR 10/6/16 -- Modify to plot all *.eig_fix files
%
% JBR 1/20/17 -- Modified to include both Spheroidal & Toroidal and
%                traction eigen functions
%

clear; close all;

mbranch = 1; % mode number
overwrite = 1; % overwrite existing eigenfunction *.mat file?


parameter_FRECHET;
periods = param.periods;
TYPE = param.TYPE;
CARDID = param.CARDID;
if ( TYPE == 'T') 
    TYPEID = param.TTYPEID;
elseif ( TYPE == 'S') 
    TYPEID = param.STYPEID;
end
mat_name = ['./eig_mats/',param.CARDID,'.',num2str(N_modes),'.',TYPE,'.mat'];

% check if eigenfunction *.mat file already exists in ./eig_mats directory
if ~exist(mat_name,'file') || overwrite
    display('Creating new eigenfunction *.mat file');
    % CHECK FOR *.eig_fix files
    
    com = ['ls ',param.TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'_1.eig_fix | cat'];
    [status eig_fils] = system(com);
    if strcmp(eig_fils(end-25:end-1),'No such file or directory')
        disp('No *.eig_fix files')
        [eigen,eig,saveopt] = load_eigen_asc(TYPE,mbranch);
        save(mat_name,'eig',saveopt);
    else
        disp('Found *.eig_fix files')
        [eigen,eig,saveopt] = load_eigen_asc(TYPE,mbranch);
        save(mat_name,'eig',saveopt);
    end
else
    display('Loading eigenfunctions from *.mat');
    load(mat_name); %loads eigenfunction structure called 'eig'
end

%% Plot eigenfunctions
fig1 = figure(1); clf;
set(gcf,'color','w','position',[121   192   451   506]);
clr = jet(length(periods));
%             clr = copper(num_ll);
            %clr = lines(num_ll);
            %clr = copper(1);
for iper = 1:length(periods)
%     [~,Iper] = min(abs([eig.nn(mbranch+1).ll(:).per]-periods(iper)));   
    %        if last_nn >= (size(eig.ll(2).nn,2) - 1)           
    r = 6371 - eig(1).r/1000;
    wl = eig(iper).wl;
    wp = eig(iper).wp;
    per = eig(iper).per;

    %subplot(1,2,1)
    %            plot(wl,r,'linewidth',2); hold on;
    plot(wl,r,'linewidth',2,'color',clr(iper,:),'linewidth',3); hold on;
    %            hold off;
    axis ij
    %            set(gca,'linewidth',2);
    title(sprintf('Eigenfunctions: nn = %d',mbranch),'fontsize',12);
    xlabel('wl');
    ylabel('depth (km)');
    %              xlim([-1 1]);
    ylim([0 400]);
    xlim([-0.1 0.1]);

    lgd{iper}=[num2str(periods(iper)),'s'];

    %            pause;
end

%subplot(1,2,1); hold on;
legend(lgd,'location','EastOutside','box','off');
plot([0,0],[0,3000],'--k','linewidth',2);
set(gca,'fontsize',16,'linewidth',2);

% subplot(1,2,2); hold on;
% legend(lgd,'location','southeast');
% plot([0,0],[0,3000],'--k');
% get(h,'interpreter');
% set(h,'interpreter','none');

CARDID = param.CARDID;
TYPEID = param.TTYPEID;
EIGPATH = param.eigpath;
%print('-painters','-dpdf','-r400',[EIGPATH,CARDID,'.',TYPEID,'.',num2str(j),'mod.',num2str(N_modes),'_fix.pdf']);
save2pdf([EIGPATH,CARDID,'.',TYPEID,'.',num2str(j),'mod.',num2str(N_modes),'_fix.pdf'],fig1,1000);

delete([param.eigpath,'*_fix.asc'],[param.eigpath,'*.asc']);