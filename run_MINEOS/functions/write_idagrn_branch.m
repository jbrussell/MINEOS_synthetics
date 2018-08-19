%% Write idagrn driver
% JBR 7/18
%
function write_idagrn_branch(TYPE,CARD,BR,EVT,STA,LENGTH_HR,DT)


parameter_FRECHET;

CARDTABLE = CARDTABLE;
CARDID = param.CARDID;

CARDPATH = param.CARDPATH;
FRECHETPATH = param.frechetpath;
TABLEPATH = param.TABLEPATH;

if strcmp(TYPE,'T') == 1
    disp('Toroidal!');
    
    RUNFILE = 'run_idagrn_branch.t';
    TYPEID = param.TTYPEID;
    
elseif strcmp(TYPE,'S') == 1
    disp('Spheroidal!');
    
    RUNFILE = 'run_idagrn_branch.s';
    TYPEID = param.STYPEID;
    
else
    disp('No TYPE recognized!');
    
end
TABLE = [CARDTABLE,CARD,'.',TYPEID,'.table'];
TABLE_MASK = [CARDTABLE,CARD,'.',TYPEID,'.',num2str(BR),'.mask'];

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
fprintf(fid,'.true.\n');
fprintf(fid,'%s\n',TABLE_MASK);
fprintf(fid,'1\n');
fclose(fid);





    