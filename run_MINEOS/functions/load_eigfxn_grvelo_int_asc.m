% load_frechcv_asc.
% loads frechet kernels from ascii files
%
% JBR 10/11/16
%

function [AMP] = load_eigfxn_grvelo_int_asc(card,TYPE)

parameter_FRECHET;

if strcmp(TYPE,'T') == 1
    disp('Toroidal!');
    
    TYPEID = param.TTYPEID;
    
elseif strcmp(TYPE,'S') == 1
    disp('Spheroidal!');
    
    TYPEID = param.STYPEID;
    
else
    disp('No TYPE recognized!');
    
end

ampfile = [param.TABLEPATH,'/',card,'/tables/',card,'.',TYPEID,'.eigfxn_grvelo'];

%     1. requested period
%     2. closest period in the mode catalog 
%     3. eigenfunction U evaluated at the surface, r=6371km
%     4. group velocity G
%     5. phase velocity c
%     6. eigenfunction V evaluated at the surface, r=6371km
fid = fopen(ampfile,'r');
C = textscan(fid,'%f%f%f%f%f%f%f');
fclose(fid);
    
if strcmp(TYPE,'S') == 1
    AMP.pers_want = C{1};
    AMP.periods = C{2};
    AMP.U_0 = C{3};
    AMP.grv = C{4};
    AMP.phv = C{5};
    AMP.V_0 = C{6};
    AMP.I_0 = C{7};
elseif strcmp(TYPE,'T') == 1
end


end