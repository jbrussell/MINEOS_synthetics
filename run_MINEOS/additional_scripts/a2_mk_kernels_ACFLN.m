%% Drive to calculate azimuthal anisotropy kernals 
% This involves calling fortran programs plot_wk, frechet_G, frechet_cvG,
% frechet_asc
% pylin.patty 11/2014
% G kernals are same as L kernals, which is written in fortran files,
% Frechet_G nad Frechet_cvG. it's for sherical modes only now. 
% pylin.patty 01/2015
%
% 2/19/2017 JBR - now handles ACFLN kernels
%
% !!! SHOULD FIRST RUN MINEOS TO GENERATE MODE TABLES !!!
% (Should have all files for every mode and branch in the */tables/
% directory in MODE_tables)
%
% 9/11/2017 JBR - No prefactor
%               - Kernels must be multiplied by c/U * rho * V^2
%
clear; close all;

parameter_FRECHET;
branch = 0; % Fundamental -> 0

is_deletefrech = 1; % Delete the .frech files to save space?

TYPE = param.TYPE;
CARDID = param.CARDID;

titlename = '';

if ( TYPE == 'T') 
    TYPEID = param.TTYPEID;
elseif ( TYPE == 'S') 
    TYPEID = param.STYPEID;
end

model_depth = 400; %80; %80 km

periods = param.periods;
FRECHETPATH = param.frechetpath;

is_frech_x = 0; % 1 => scale ax; 0 => autoscale
ylims = [0 100]; %[0 400];
%ylims = [0 700];
xlims = [0 1e-17]*1000*4;



isfigure = 1;

%% Set path to executables
setpath_plotwk;

%% Change environment variables to deal with gfortran
setenv('GFORTRAN_STDIN_UNIT', '5') 
setenv('GFORTRAN_STDOUT_UNIT', '6') 
setenv('GFORTRAN_STDERR_UNIT', '0')


%% run plot_wk on the table_hdr file to generate the branch file
write_plotwk(TYPE,CARDID);

com = ['cat run_plotwk.',lower(TYPE),' | plot_wk'];
[status,log] = system(com);

if status ~= 0     
    error( 'something is wrong at plot_wk')
end

%% run "frechet" to generate the frechet file 
NDISC = 0;
ZDISC = [];
if ( TYPE == 'T') 
    TYPEID = param.TTYPEID;
elseif ( TYPE == 'S') 
    TYPEID = param.STYPEID;
end

com = ['ls ',param.TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'_1.eig_fix | cat'];
[status eig_fils] = system(com);
if strcmp(eig_fils(end-25:end-1),'No such file or directory')
    disp('Found no *.eig_fix files')
    write_frechet(TYPE,CARDID,NDISC,ZDISC)
else
    disp('Found *.eig_fix files')
    write_frech_chk(NDISC)
end
disp('Be patient! This will take ~25 s');
tic
com = ['cat run_frechet.',lower(TYPE),' | frechet_ACFLN_love > frechet.LOG'];
[status,log] = system(com);
if status ~= 0     
    error( 'something is wrong at frechet_ACFLN_love')
end
toc

%% run "frechet_cv" to generate the cv kernels
% Convert frechet to ascii
% Make CV Frechet Kernels
disp('--- Make CV Frechet Kernels ---');


write_frechcv(TYPE,CARDID,branch)

com = ['cat run_frechcv.',lower(TYPE),' | frechet_cvG'];
[status,log] = system(com);
if status ~= 0     
    error( 'something is wrong at frechet_cvG')
end
%% Convert CV Frechet kernels to ascii with phase-velocity sensitivity
% Will do this for all periods of interest
% Set inside the setparam_MINE.m
disp('--- Convert Frechet CV to ascii ---');

    % Program writes run file for draw_frechet_gv, runs it, and reads in
    % sensitivity kernels for all periods of interest

FRECH = frechACFLN_asc(TYPE,CARDID,branch);  

%% Plot Kernels

fig1 = figure(1); clf;
%set(gcf,'position',[58   255   548   450]);
% set(gcf,'position',[58   255   916   450]);

clr = jet(length(periods));
rad = (FRECH(1).rad(end)-FRECH(1).rad)/1000;       
dr = gradient(rad);        
for iper = 1:length(FRECH) %length(periods)
    lgd{iper} = [num2str(FRECH(iper).per),' s'];

    if ( TYPE == 'T') 
        set(gcf,'position',[58   255   916   450],'color','w');

        % L
        subplot(1,3,1);
        set(gca,'linewidth',2);
        h(iper) = plot(FRECH(iper).L .* dr,rad,'linewidth',2,'color',clr(iper,:)); hold on;
        axis ij
        title(['L (2\theta)'],'fontsize',15);
        ylabel('Depth (km)','fontsize',15);
        ylim(ylims);
        %xlim([0 2e-5]);
        if  is_frech_x
        	lim(xlims*30);
        end
        set(gca,'fontsize',15)

        % N
        subplot(1,3,2);
        set(gca,'linewidth',2);
        plot(FRECH(iper).N .* dr,rad,'linewidth',2,'color',clr(iper,:)); hold on;
        axis ij
        title(['N (4\theta)'],'fontsize',15);
        ylabel('Depth (km)','fontsize',15);
        ylim(ylims);
%         xlim([1e-16 1e-11]);
        if  is_frech_x
            xlim(xlims*200);
        end
        set(gca,'fontsize',15)

        % RHO
        subplot(1,3,3);
        set(gca,'linewidth',2);
        plot(FRECH(iper).rho .* dr,rad,'linewidth',2,'color',clr(iper,:)); hold on;
        axis ij
        title(['\rho Kernels'],'fontsize',15);
        ylabel('Depth (km)','fontsize',15);
        ylim(ylims);
%         xlim([1e-16 1e-11]);
        if  is_frech_x
            xlim(xlims);
        end
        set(gca,'fontsize',15)



    elseif ( TYPE == 'S') 
        set(gcf,'position',[10         234        1174         471]);

        % A
        subplot(1,5,1);
        set(gca,'linewidth',2);
        h(iper) = plot(FRECH(iper).A .* dr,rad,'linewidth',2,'color',clr(iper,:)); hold on;
        axis ij
        title(['A (2\theta & 4\theta)'],'fontsize',15);
        ylabel('Depth (km)','fontsize',15);
        ylim(ylims);
%         xlim([1e-16 1e-11]);
        if  is_frech_x
            xlim(xlims*15);
        end
        set(gca,'fontsize',15)

        % C
        subplot(1,5,2);
        set(gca,'linewidth',2);
        plot(FRECH(iper).C .* dr,rad,'linewidth',2,'color',clr(iper,:)); hold on;
        axis ij
        title(['C'],'fontsize',15);
        ylabel('Depth (km)','fontsize',15);
        ylim(ylims);
%         xlim([1e-16 1e-11]);
        if  is_frech_x
            xlim(xlims*5);
        end
        set(gca,'fontsize',15)

        % F
        subplot(1,5,3);
        set(gca,'linewidth',2);
        if iper == 1
            plot([0 0],ylims,'--k'); hold on;
        end
        plot(FRECH(iper).F .* dr,rad,'linewidth',2,'color',clr(iper,:)); hold on;
        axis ij
        title(['F (2\theta)'],'fontsize',15);
        ylabel('Depth (km)','fontsize',15);
        ylim(ylims);
%         xlim([1e-16 1e-11]);
        if  is_frech_x
            xlim(xlims*5);
        end
        set(gca,'fontsize',15)

        % L
        subplot(1,5,4);
        set(gca,'linewidth',2);
        plot(FRECH(iper).L .* dr,rad,'linewidth',2,'color',clr(iper,:)); hold on;
        axis ij
        title(['L (2\theta)'],'fontsize',15);
        ylabel('Depth (km)','fontsize',15);
        ylim(ylims);
%         xlim([1e-16 1e-11]);
        if  is_frech_x
            xlim(xlims*20);
        end
        set(gca,'fontsize',15)

        % RHO
        subplot(1,5,5);
        set(gca,'linewidth',2);
        plot(FRECH(iper).rho .* dr,rad,'linewidth',2,'color',clr(iper,:)); hold on;
        axis ij
        title(['\rho'],'fontsize',15);
        ylabel('Depth (km)','fontsize',15);
        ylim(ylims);
%         xlim([1e-16 1e-11]);
        if  is_frech_x
            xlim(xlims);
        end
        set(gca,'fontsize',15)
    end
    %display(FRECH(iper).per)
    %pause;
end

if ( TYPE == 'T')
    subplot(1,3,1); hold on;
    legend(h,lgd,'location','southeast','fontsize',15);
elseif ( TYPE == 'S') 
    subplot(1,5,1); hold on;
    legend(h,lgd,'location','southeast','fontsize',15);
    % plot([0,0],[0,3000],'--k');
end

FRECHETPATH = param.frechetpath;
delete(['run_plotwk.',lower(TYPE)],['run_frechcv.',lower(TYPE)],['run_frechet.',lower(TYPE)],['run_frechcv_asc.',lower(TYPE)]);
save2pdf([FRECHETPATH,'CARD_ACFLN_kernels_',lower(TYPE),'_',CARDID,'_b',num2str(branch),'.',num2str(N_modes),'_',num2str(periods(1)),'_',num2str(periods(end)),'s.pdf'],fig1,1000)


%    savefile = [CARD,'_fcv.mat'];
%    save(savefile,'FRECH_T','FRECH_S');
% Change the environment variables back to the way they were
setenv('GFORTRAN_STDIN_UNIT', '-1') 
setenv('GFORTRAN_STDOUT_UNIT', '-1') 
setenv('GFORTRAN_STDERR_UNIT', '-1')

delete('*.LOG');
if is_deletefrech
    delete([param.frechetpath,'*.fcv.*'],[param.frechetpath,'*.frech'])
end
