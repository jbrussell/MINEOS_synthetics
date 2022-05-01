% write_MINE_qmod
function [] = write_MINEOS_qmod(rad_km,qmu,qkap,CARD)

% CARD : Name of output card file
% ncard : input card structure to write out

%% Write qmod

fid=fopen(CARD,'w');
fprintf(fid, '%13i%10s%10s\n',length(qmu),'shear','bulk');
for id = 1:length(qmu)
    fprintf(fid,'%13.2f%10.2f%10.2f\n'...
        ,[rad_km(id) qmu(id) qkap(id)]);
end
fclose(fid);

end
