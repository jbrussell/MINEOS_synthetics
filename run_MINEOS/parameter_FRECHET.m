%% Setup Parameters for running MINEOS to calculate senstivity kernels, dispersion, and synthetics
% NJA, 2014
% pylin.patty, 2014
%
%9/14/15: Josh Russell
%         Added path for eigenfunctions 'param.eigpath'
%10/8/15: Josh Russell
%         Added path for dispersion curves 'param.disperspath'
%7/21/18: JBR
%         Can now calculate synthetics using idagrn
%

%clear all;
addpath('./functions');

path2runMINEOS = './';

% Mineos table parameters
maxN = 400000; % Estimate of max number of modes
minF = 0;
maxF = 200.05; %200.05; %10.1; %200.05; %250.05; %333.4; %500.05; %200.05; %%150.05; %50.05;
minL = 0;
maxL = 50000; %50000; %6000;
N_modes = 2; % <0 uses all mode branches, 1=fundamental only -------- JOSH 8/22/15
param.CARDID = 'pa5_5km'; %'Nomelt_taper_aniso_constxicrman_etaPREM_constxilays'; %'pa5_5km';

ch_mode = 0; % mode branch to check for missed eigenfrequencies 0 => T0 ------- JOSH 10/7/15
%param.TYPE = 'T';

% (1 => yes, 0 => no)
SONLY = 0; %Spheroidal modes? (RAYLEIGH)
TONLY = 1; %Toroidal modes? (LOVE)

%% Parameters for idagrn synthetics
LENGTH_HR = 1.0; %1.0; % length of seismogram in hours
DT = 1.0; % 1/samplerate
eventfile = 'evt_0000001';
stationfile = 'stations.stn';


% Setup idagrn paths
param.IDAGRN = [path2runMINEOS,'/IDAGRN/'];
param.EVTPATH = [param.IDAGRN,'EVT_FILES/',eventfile];
param.STAPATH = [param.IDAGRN,'STATION/',stationfile];
param.SYNTH_OUT = [param.IDAGRN,'SYNTH/',param.CARDID,'_b',num2str(N_modes),'/',eventfile,'/'];
if ~exist(param.SYNTH_OUT)
    mkdir(param.SYNTH_OUT);
end

%%
if SONLY == 1 && TONLY == 0
    param.TYPE = 'S';
elseif SONLY == 0 && TONLY == 1
    param.TYPE = 'T';
else
    error('Choose SONLY or TONLY, not both');
    
end

% Setup Parameters for Initial Model
param.CARD = [param.CARDID,'.card'];
param.CARDPATH  = [path2runMINEOS,'/CARDS/'];
param.TABLEPATH = [path2runMINEOS,'/MODE/TABLES/'];
param.MODEPATH  = [path2runMINEOS,'/MODE/TABLES/MODE.in/'];
param.RUNPATH = pwd;

%% create dir for output MINEOS automatically, doesn't need to be changed.
CARDTABLE = [param.TABLEPATH,param.CARDID,'/tables/'];
if ~exist(CARDTABLE)
    mkdir([param.TABLEPATH,param.CARDID])
    mkdir(CARDTABLE)
end

%% setup Parameters for kernals
param.frechetpath = [path2runMINEOS,'/MODE/FRECHET/',param.CARDID,'/'];

if ~exist(param.frechetpath) 
    mkdir(param.frechetpath)
end

%% setup Parameters for eigenfunctions
param.eigpath = [path2runMINEOS,'/MODE/EIGEN/',param.CARDID,'/'];

if ~exist(param.eigpath) 
    mkdir(param.eigpath)
end

%% setup Parameters for Dispersion
param.disperspath = [path2runMINEOS,'/MODE/DISPERSION/',param.CARDID,'/'];

if ~exist(param.disperspath) 
    mkdir(param.disperspath)
end

%% Turn on if only want to calculate S or T or both for mineous
param.SMODEIN = ['s.mode',num2str(floor(minF)),'_',num2str(floor(maxF)),'_b',num2str(N_modes)];
param.STYPEID = ['s',num2str(floor(minF)),'to',num2str(floor(maxF))];
param.TMODEIN = ['t.mode',num2str(floor(minF)),'_',num2str(floor(maxF)),'_b',num2str(N_modes)];
param.TTYPEID = ['t',num2str(floor(minF)),'to',num2str(floor(maxF))];%'t0to150';

%% plotting kernels
% param.periods = [20 25 32 40 50 60 80 100 115 130 145];
%param.periods = [10 12 14 16 18 20 25 32 40 50 60 80 100 115 130 145];
%param.periods = [5 7 10 12 14 16 18 20 25 32 40 50 60 80 100 115 130 145];
%param.periods = [3 5 7 10 12 14 16 18 20 25 32 40 50 60 80 100 115 130 145];
%param.periods = [3 5 7 10];
%param.periods = [4 5 6 7 8 9 10 12 14 16 18 20 25 32 40 50 60 80 100]; %[4 5 6 7 8 9 10]; %[4 5 7 10 12 14 16 18 20 25 32 40 50 60 80 100]; %[4 5 6 7 8 9 10]; %[2 4 6 8 10 12]; %[20 25 32 40 50 60 80 100];
%param.periods = [ 12 ]

% param.periods = flip([7.5000    6.6667    6.0000    5.4545    5.0000]);

%param.periods = flip([30.0000   25.7143   22.5000   20.0000   18.0000   16.3636  15.0000   13.8462   12.8571   12.0000]);
% param.periods = flip([30.0000   25.7143   22.5000   20.0000   18.0000   16.3636  15.0000   13.8462   12.8571   12.0000 7.5000    6.6667    6.0000    5.4545    5.0000]);
 %param.periods = flip([30.0000   25.7143   22.5000   20.0000   18.0000   16.3636  15.0000   13.8462   12.8571   12.0000    10.0000    8.5714   7.5000    6.6667    6.0000    5.4545    5.0000]);
%param.periods = flip([10.0000    8.5714   7.5000    6.6667    6.0000    5.4545    5.0000]);

%param.periods = [82, 212]; % Rayleigh

% % for Helen
% param.periods = [40 50 60 80 100 115 130 145 160];

% % for Zach
% param.periods = round(logspace(log10(10),log10(200),12));

param.periods = round(logspace(log10(5),log10(200),15));
