%% Write idagrn driver
% JBR 7/18
%
function write_idagrn(TYPE,CARD,EVT,STA,LENGTH_HR,DT)


parameter_FRECHET;

CARDTABLE = CARDTABLE;
CARDID = param.CARDID;

CARDPATH = param.CARDPATH;
FRECHETPATH = param.frechetpath;
TABLEPATH = param.TABLEPATH;

if strcmp(TYPE,'T') == 1
    disp('Toroidal!');
    
    RUNFILE = 'run_idagrn.t';
    TYPEID = param.TTYPEID;
    
elseif strcmp(TYPE,'S') == 1
    disp('Spheroidal!');
    
    RUNFILE = 'run_idagrn.s';
    TYPEID = param.STYPEID;
    
else
    disp('No TYPE recognized!');
    
end
TABLE = [CARDTABLE,CARD,'.',TYPEID,'.table'];

%% Get event information
com = sprintf('awk ''{print $0}'' %s',EVT);
[log, EVTSTR] = system(com);

%% Get station information
com = sprintf('awk ''{print $0}'' %s',STA);
[log, STASTR] = system(com);

%% Write file

if exist(RUNFILE,'file') == 2
    disp('File exists! Removing it now')
    com = ['rm -f',RUNFILE];
    [status,log] = system(com);
end

fid = fopen(RUNFILE,'w');
fprintf(fid,'%s\n',TABLE);
fprintf(fid,'%.4f %.1f\n',LENGTH_HR,DT);
fprintf(fid,'35.\n');
fprintf(fid,'%s\n',CARDID);
fprintf(fid,'%s\n',TYPE);
fprintf(fid,'.true.\n');
fprintf(fid,'.false.\n');
fprintf(fid,'sa_junk\n');
fprintf(fid,'frechet_file_junk\n');
fprintf(fid,'%s',EVTSTR);
fprintf(fid,'%s',STASTR);
fclose(fid);





    