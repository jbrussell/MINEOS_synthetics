% frechcv_asc.
% Program converts binary frechet kernels in phase velocity to ascii
% versions and saves the output in a mat file
% Involves writing run files for draw_frechet_gv and running the fortran
% program
%
% JBR 11/22/2021 - read in Q kernels K_Qmu and K_Qkappa. They are
% premultiplied by the relevant term according to Dziewonski and Anderson
% 1981
%    - K_Qmu = mu * MM
%    - K_Qkappa = kappa * KK
%
%

function [FRECH] = frechQ_asc(TYPE,CARD,BRANCH)

% TYPE = 'T';
% CARD = param.CARDID;
% BRANCH = 0;

fid_log = fopen('draw_frechet_gv.LOG','w'); %JBR -- write log file

% Get useful info from parameter file
parameter_FRECHET;
CARDPATH = param.CARDPATH;
FRECHETPATH = param.frechetpath;
TABLEPATH = param.TABLEPATH;
periods = param.periods;

if strcmp(TYPE,'T') == 1
    disp('Toroidal!');
    
    RUNFILE = 'run_frechcv_asc.t';
    TYPEID = param.TTYPEID;
    
elseif strcmp(TYPE,'S') == 1
    disp('Spheroidal!');
    
    RUNFILE = 'run_frechcv_asc.s';
    TYPEID = param.STYPEID;
    
else
    disp('No TYPE recognized!');
    
end

BRID = [num2str(BRANCH)];
% if BRANCH == 0
%     BRID = '0st';
% elseif BRANCH == 1
%     BRID = '1st';
% elseif BRANCH == 2
%     BRID = '2nd';
% elseif BRANCH == 3
%     BRID = '3rd';
% else
%     disp('Branch has no name! Change it in the script')
% end

FRECHCV = [FRECHETPATH,CARD,'.',TYPEID,'.fcv.',BRID];




for ip = 1:length(periods)
    
    FRECHASC = [FRECHETPATH,CARD,'.',TYPEID,'.',BRID,'.',num2str(periods(ip))];
    
    % Write runfile for draw_frechet_gv
    fid = fopen(RUNFILE,'w');
    fprintf(fid,'%s\n',FRECHCV); %input binary file
    fprintf(fid,'%s\n',FRECHASC); %output ascii file
    fprintf(fid,'%i\n',periods(ip));
    fclose(fid);
    
    % Run draw_frechet_gv
    disp(sprintf('--- Period : %s',num2str(periods(ip))));
    
    if exist(FRECHASC,'file') == 2
    %disp('File exists! Removing it now')
    com = ['rm -f ',FRECHASC];
    [status,log] = system(com);
    end
    
    com = sprintf('cat %s | draw_frechet_gv',RUNFILE);
    [status,log] = system(com);
    fprintf(fid_log,log); % JBR
%    log
    dT = abs(periods(ip)-str2num(log(end-10:end)));
    if ( dT < 1.5)
        disp(sprintf('Find closest period: %s',log(end-10:end)));
    else
        disp('the closest period is FAR AWAY from period of interest')
        disp(sprintf('Find closest period: %s',log(end-10:end)));
        
    end
    
    % Load in frechet files for each period
    
    % spheroidal, no aniso: 1=Vs,2=Vp,3=rho
    % spheroidal, aniso: 1=Vsv,2=Vpv,3=Vsh,4=Vph,5=eta,6=rho
    % toroidal, no aniso: 1=Vs,2=rho
    % toroidal, aniso: 1=Vsv,2=Vsh,3=rho
    % disp(FRECHASC)
    fid = fopen(FRECHASC,'r');
    C = textscan(fid,'%f%f%f%f%f%f%f');
    fclose(fid);
    
    if strcmp(TYPE,'S') == 1
        FRECH(ip).per = periods(ip);
        FRECH(ip).rad = C{1};
        FRECH(ip).K_qmu = C{2};
        FRECH(ip).K_qkappa = C{3};
        FRECH(1).qmu = C{4}; % qmu from qmod interpolated to card knots
        FRECH(1).qkappa = C{5}; % qkappa from qmod interpolated to card knots
%         FRECH(ip).F = C{6}; % dummy
        FRECH(ip).K_rho = C{7};
    elseif strcmp(TYPE,'T') == 1
        FRECH(ip).per = periods(ip);
        FRECH(ip).rad = C{1};
        FRECH(ip).K_qmu = C{2};
        FRECH(1).qmu = C{3}; % qmu from qmod interpolated to card knots
        FRECH(ip).K_rho = C{4};
    end
end

fclose(fid_log); % JBR
