%% Driver to Calculate Frechet Kernels in terms of Phase Velocity
% NJA, 2014
% 
% This involves calling fortran programs plot_wk, frechet, frechet_gv,
% frechet_pv
% pylin.patty 2014
%
% JOSH 8/25/2015 --- only plot T and plot CARD file
%


clear; close all;

parameter_FRECHET;
branch = 1; % Fundamental -> 0

TYPE = param.TYPE;
CARDID = param.CARDID;

titlename = 'Patty Isotropic';

if ( TYPE == 'T') 
    TYPEID = param.TTYPEID;
elseif ( TYPE == 'S') 
    TYPEID = param.STYPEID;
end

periods = param.periods;

yaxis = [0 350]; %[0 100]; %[0 350];

is_frech_x = 1; % 1 => scale ax; 0 => autoscale
frech_x = [0 2e-8]; %[0 2.0e-7]; %[0 2.0e-8];
%frech_x = [-2.5e-8 2.5e-8];

isfigure = 1;

card_comp = 'pa5_5km'; %'Nomelt3inttaper'; %'prem_oceanint'; %'iasp91_jim'; %'pa5iso'; %'pa5'; %'iasp91';

%% Set path to executables
setpath_plotwk;

%% Change environment variables to deal with gfortran
setenv('GFORTRAN_STDIN_UNIT', '5') 
setenv('GFORTRAN_STDOUT_UNIT', '6') 
setenv('GFORTRAN_STDERR_UNIT', '0')


%% run plot_wk on the table_hdr file to generate the branch file
write_plotwk(TYPE,CARDID);

com = ['cat run_plotwk.',lower(TYPE),' | plot_wk > plot_wk.LOG'];
[status,log] = system(com);

if status ~= 0     
    error( 'something is wrong at plot_wk')
end


%% run "frechet" to generate the frechet file 
NDISC = 0;
ZDISC = [];

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
com = ['cat run_frechet.',lower(TYPE),' | frechet > frechet.LOG'];
[status,log] = system(com);
if status ~= 0     
    error( 'something is wrong at frechet')
end
toc

%% run "frechet_cv" to generate the cv kernels
% Convert frechet to ascii
% Make CV Frechet Kernels
disp('--- Make CV Frechet Kernels ---');


write_frechcv(TYPE,CARDID,branch)

com = ['cat run_frechcv.',lower(TYPE),' | frechet_cv > frechet_cv.LOG'];
[status,log] = system(com);
if status ~= 0     
    error( 'something is wrong at frechet_cv')
end

%% load CARD file (vmod)

card_comp_PATH = [param.CARDPATH,card_comp,'.card']; %--JBR
fid = fopen(card_comp_PATH);

%skip 3 line header
for i=1:3
    fgetl(fid);
end

%read remaining lines of file and build radius, vp, and vs vectors
pa5file = [];
while ~feof(fid)
    line = fgetl(fid);
    line = line(1:80);
    pa5file = [pa5file; line];
end
fclose(fid);

r_0core_pa5 = str2num(pa5file(:,1:7))/1000;
r_pa5 = 6371 - r_0core_pa5;

vp_pa5 = str2num(pa5file(:,19:26))/1000;
%vs_pa5 = str2num(pa5file(:,29:35))/1000; %vsv
vs_pa5 = str2num(pa5file(:,65:71))/1000; %vsh
%%%%%%% now load other card
CARD = param.CARD;
CARDPATH = param.CARDPATH;
FULLPATH = [CARDPATH,CARD];

fid = fopen(FULLPATH);

%skip 3 line header
for i=1:3
    fgetl(fid);
end

%read remaining lines of file and build radius, vp, and vs vectors
% cardfile = [];
% while ~feof(fid)
%     line = fgetl(fid);
%     line = line(1:80);
%     cardfile = [cardfile; line];
% end
% fclose(fid);
% 
% r_0core = str2num(cardfile(:,1:7))/1000;
% r = 6371 - r_0core;
% 
% vp = str2num(cardfile(:,19:26))/1000; %vpv
% %vs = str2num(cardfile(:,29:35))/1000; %vsv
% vs = str2num(cardfile(:,65:71))/1000; %vsh
ncard = textscan(fid, '%f%f%f%f%f%f%f%f%f');
fclose(fid);

ncard_temp = ncard;
R = ncard{1};
RHO = ncard{2};
VPV = ncard{3};
VSV = ncard{4};
QKAPPA = ncard{5};
QSHEAR = ncard{6};
VPH = ncard{7};
VSH = ncard{8};
eta = ncard{9};

vs = VSV/1000;
vp = VPV/1000;
r = 6371-R/1000;
%% Convert CV Frechet kernels to ascii with phase-velocity sensitivity
% Will do this for all periods of interest

disp('--- Convert Frechet CV to ascii ---');

    % Program writes run file for draw_frechet_gv, runs it, and reads in
    % sensitivity kernels for all periods of interest
    
if ( TYPE == 'S') 
    FRECH_S = frechcv_asc(TYPE,CARDID,branch);
    if isfigure
        fig1 = figure(62); set(gcf, 'Color', 'w');
        clf
         CC=lines(length(periods));
         %CC=copper(length(periods));
            subplot(1,2,1)
            %subplot(1,3,2)
            plot(vs_pa5,r_pa5,'--','linewidth',2,'color',[.5 .5 .5]); hold on;
            plot(vs,r,'linewidth',2,'color',[0 0 1]);
            set(gca,'Ydir','reverse');
            ylim([yaxis]);
            %ylim([0 100]);
            %xlim([3 5]);
            xlim([0 5.2]);
            title('Vs');
            leg = legend(card_comp,param.CARDID,'location','southwest');
            get(leg,'interpreter');
            set(leg,'interpreter','none');

    %         subplot(1,3,1)
    %         plot(vs_pa5,r_pa5,'--','linewidth',2,'color',[.5 .5 .5]); hold on;
    %         plot(vs,r,'linewidth',2,'color',[0 0 1]);
    %         set(gca,'Ydir','reverse');
    %         ylim([4 15]);
    %         %ylim([0 100]);
    %         xlim([0 5]);
    %         title('Vs');

            for ip = 1:length(periods)
                subplot(1,2,2)
                %subplot(1,3,3)
                hold on
                %if periods(ip) > 18
                plot(FRECH_S(ip).vsv,(6371000-FRECH_S(ip).rad)./1000,'-k','linewidth',2,'color',CC(ip,:))
                %end
                title('FRECH S - Vsv','fontname','Times New Roman','fontsize',12);
                lgd{ip}=[num2str(periods(ip)),' S'];
                set(gca,'YDir','reverse')
                ylim(yaxis)
                if is_frech_x == 1
                    xlim(frech_x);
                end

                box on;
    %            pause;

            end
            legend(lgd,'Location','southeast');
            %xlim([0 2e-8])



    end
    FRECHETPATH = param.frechetpath;
    delete(['run_plotwk.',lower(TYPE)],['run_frechcv.',lower(TYPE)],['run_frechet.',lower(TYPE)],['run_frechcv_asc.',lower(TYPE)]);
    %print('-painters','-dpdf','-r400',[FRECHETPATH,'CARD_Vs_kernels_',lower(TYPE),'_',CARDID,'_b',num2str(branch),'.',num2str(N_modes),'.pdf']);
    export_fig(fig1,[FRECHETPATH,'CARD_Vsv_kernels_',lower(TYPE),'_',CARDID,'_b',num2str(branch),'.',num2str(N_modes),'.pdf'],'-pdf','-painters');
    %save2pdf([FRECHETPATH,'Vs_kernels_',lower(TYPE),'_',CARDID,'.pdf'],62,300)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ( TYPE == 'T')
    FRECH_T = frechcv_asc(TYPE,CARDID,branch);  

    if isfigure
        fig1 = figure(62); set(gcf, 'Color', 'w'); 
%         set(gcf,'position',[360   165   560   532]);
        set(gcf,'position',[112   169   830   532]);
        clf
%          CC=jet(length(periods));
         CC=lines(length(periods));
         %CC=copper(length(periods));
            ax1 = subplot(1,3,1);
            %subplot(1,3,2)
            %plot(vs_pa5,r_pa5,'--','linewidth',2,'color',[.5 .5 .5]); hold on;
%             plot(vs_pa5,r_pa5,'linewidth',2,'color',[1 0 0]); hold on;
            plot(vs,r,'linewidth',3,'color',[0 0 0]);
            set(gca,'Ydir','reverse','linewidth',2,'YMinorTick','on','XMinorTick','on');
            ylim([yaxis]);
            %ylim([0 100]);
            %xlim([3 5]);
            xlim([3 5.2]);
%             title('Velocity Model','fontsize',18);
            xlabel('V_{SH} (km/s)','fontsize',18);
            ylabel('Depth (km)','fontsize',18);
            %leg = legend(card_comp,param.CARDID,'location','southwest');
            %leg = legend('Prem Ocean','Pa5','location','southwest');
            %get(leg,'interpreter');
            %set(leg,'interpreter','none');
            set(gca,'fontsize',18);
            dx = 0.06;
            ax1.Position = [ax1.Position(1:2) ax1.Position(3)+dx ax1.Position(4)];
            

    %         subplot(1,3,1)
    %         plot(vs_pa5,r_pa5,'--','linewidth',2,'color',[.5 .5 .5]); hold on;
    %         plot(vs,r,'linewidth',2,'color',[0 0 1]);
    %         set(gca,'Ydir','reverse');
    %         ylim([4 15]);
    %         %ylim([0 100]);
    %         xlim([0 5]);
    %         title('Vs');

            for ip = 1:length(periods)
                ax2 = subplot(1,3,2);
                %subplot(1,3,3)
                axis tight;
                hold on
                %if periods(ip) > 18
                plot(FRECH_T(ip).vsh,(6371000-FRECH_T(ip).rad)./1000,'-k','linewidth',3,'color',CC(ip,:))
                %end
%                 title(titlename,'fontsize',18);
                lgd{ip}=[num2str(periods(ip)),'s'];
                set(gca,'YDir','reverse','linewidth',2,'YMinorTick','on','XMinorTick','on')
                ylim(yaxis)
                xlabel('dc/dV_{SH}');
                if is_frech_x == 1
                    xlim(frech_x);
                end

                box on;
    %            pause;

            end

            %legend(lgd,'Location','eastoutside');
            %xlim([0 2e-8])
            fig=gcf;
            set(gca,'fontsize',18);
            ax2.Position = [ax1.Position(1)+ax1.Position(3)+0.1 ax2.Position(2) ax2.Position(3)+dx ax2.Position(4)];
            legend(lgd,'position',[ax2.Position(1)+ax2.Position(3)+0.08 0.5 0 0],'box','off');


    end
    FRECHETPATH = param.frechetpath;
    delete(['run_plotwk.',lower(TYPE)],['run_frechcv.',lower(TYPE)],['run_frechet.',lower(TYPE)],['run_frechcv_asc.',lower(TYPE)]);
    %print('-painters','-dpdf','-r400',[FRECHETPATH,'CARD_Vs_kernels_',lower(TYPE),'_',CARDID,'_b',num2str(branch),'.',num2str(N_modes),'.pdf']);
    %export_fig(fig1,[FRECHETPATH,'CARD_Vs_kernels_',lower(TYPE),'_',CARDID,'_b',num2str(branch),'.',num2str(N_modes),'.pdf'],'-pdf','-painters');
    save2pdf([FRECHETPATH,'CARD_Vs_kernels_',lower(TYPE),'_',CARDID,'_b',num2str(branch),'.',num2str(N_modes),'_',num2str(periods(1)),'_',num2str(periods(end)),'s.pdf'],fig1,1000)
end
%    savefile = [CARD,'_fcv.mat'];
%    save(savefile,'FRECH_T','FRECH_S');
% Change the environment variables back to the way they were
setenv('GFORTRAN_STDIN_UNIT', '-1') 
setenv('GFORTRAN_STDOUT_UNIT', '-1') 
setenv('GFORTRAN_STDERR_UNIT', '-1')

delete('*.LOG');
delete([param.frechetpath,'*.fcv.*'],[param.frechetpath,'*.frech'])