%% Write frechet_gv driver
% NJA, 2014
% Must designate mode branch of interest (0 = fundamental)
% make TYPEID as a parameter in parameter_FRECHET
% pylin.patty 2015/01
function write_frechgv(TYPE,CARD,BR)


parameter_FRECHET;

CARDPATH = param.CARDPATH;
% FRECHETPATH = param.frechetpath;
FRECHETPATH = [param.frechet,CARD,'/'];
TABLEPATH = param.TABLEPATH;

if strcmp(TYPE,'T') == 1
    disp('Toroidal!');
    
    RUNFILE = 'run_frechgv.t';
    TYPEID = param.TTYPEID;
    
elseif strcmp(TYPE,'S') == 1
    disp('Spheroidal!');
    
    RUNFILE = 'run_frechgv.s';
    TYPEID = param.STYPEID;
    
else
    disp('No TYPE recognized!');
    
end

BRID = [num2str(BR)];
% if BR == 0
%     BRID = '0st';
% elseif BR == 1
%     BRID = '1st';
% elseif BR == 2
%     BRID = '2nd';
% elseif BR == 3
%     BRID = '3rd';
% else
%     disp('Branch has no name! Change it in the script')
% end

% QMOD = [CARDPATH,CARD,'.qmod'];
% BRANCH = [DATAPATH,CARD,'.',TYPEID,'.table_hdr.branch'];
FRECH = [FRECHETPATH,CARD,'.',TYPEID,'.frech'];
FRECHGV = [FRECHETPATH,CARD,'.',TYPEID,'.fgv.',BRID];

if exist(FRECHGV,'file') == 2
    disp('File exists! Removing it now')
    com = ['rm -f',FRECHGV];
    [status,log] = system(com);
end

fid = fopen(RUNFILE,'w');
fprintf(fid,'%s\n',FRECH);
fprintf(fid,'%i\n',BR);
fprintf(fid,'%s\n',FRECHGV);
fclose(fid);

    