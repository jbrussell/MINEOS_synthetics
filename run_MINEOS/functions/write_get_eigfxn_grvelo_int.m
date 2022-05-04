%% Write run file for frechet.f
% NJA, 2014
% make TYPEID as a parameter in parameter_FRECHET
% pylin.patty 2015/01
%
function write_get_eigfxn_grvelo_int(periods,ind_surf)

parameter_FRECHET;
CARDPATH  = param.CARDPATH;
TABLEPATH = param.TABLEPATH;
% Set number of discontinuities to add

parameter_FRECHET;
TYPE = param.TYPE;

if strcmp(TYPE,'S') == 1
    TYPEID = param.STYPEID;
    RUNFILE_eigfxn_grvelo = 'run_get_eigfxn_grvelo.s';
elseif strcmp(TYPE,'T') == 1
    TYPEID = param.TTYPEID;
    RUNFILE_eigfxn_grvelo = 'run_get_eigfxn_grvelo.t';
end

EIG = [CARDTABLE,param.CARDID,'.',TYPEID,'.eig'];
AMP = [CARDTABLE,param.CARDID,'.',TYPEID,'.eigfxn_grvelo'];

% Write eig_recover driver
fid=fopen(RUNFILE_eigfxn_grvelo,'w');
fprintf(fid,'%s\n',EIG);
fprintf(fid,'%d\n',length(periods));
for ip = 1:length(periods)
    fprintf(fid,'%f\n',periods(ip));
end
fprintf(fid,'%s\n',AMP);
fprintf(fid,'%d\n',ind_surf);
fclose(fid);
