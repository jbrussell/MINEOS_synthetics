function ev = read_cmtsolution(fname)
% ev = read_cmtsolution(fname)
% Read Global CMT solutions in CMTSOLUTION format (blocks separated by
% blank lines; lines starting with # are comments).
% ev(k): id, region, lat, lon, dep (km, centroid), M = [Mrr Mtt Mpp Mrt
% Mrp Mtp] (dyne-cm), Mw, hdur (s), tshift (s), and the idagrn6 inputs
% fac and mt (mt in units of 1e20 dyne-cm / fac).
txt = fileread(fname);
lines = regexp(txt, '\r?\n', 'split');
ev = struct('id', {}, 'region', {}, 'lat', {}, 'lon', {}, 'dep', {}, 'M', {}, ...
            'Mw', {}, 'hdur', {}, 'tshift', {}, 'fac', {}, 'mt', {});
k = 0; cur = [];
for i = 1:numel(lines)
    l = lines{i};
    if isempty(strtrim(l)) || l(1) == '#', continue; end
    if isempty(strfind(l, ':'))                 % hypocenter line starts a block
        k = k + 1; cur = struct();
        w = strsplit(strtrim(l));
        cur.region = strjoin(w(12:end), ' ');
        ev(k).region = cur.region;
        continue
    end
    [key, val] = strtok(l, ':'); key = strtrim(key); val = strtrim(val(2:end));
    switch key
        case 'event name', ev(k).id = val;
        case 'time shift', ev(k).tshift = str2double(val);
        case 'half duration', ev(k).hdur = str2double(val);
        case 'latitude', ev(k).lat = str2double(val);
        case 'longitude', ev(k).lon = str2double(val);
        case 'depth', ev(k).dep = str2double(val);
        case {'Mrr','Mtt','Mpp','Mrt','Mrp','Mtp'}
            j = find(strcmp(key, {'Mrr','Mtt','Mpp','Mrt','Mrp','Mtp'}));
            ev(k).M(j) = str2double(val);
    end
end
for k = 1:numel(ev)
    M = ev(k).M;
    M0 = sqrt(0.5*(sum(M(1:3).^2) + 2*sum(M(4:6).^2)));
    ev(k).Mw = 2/3*log10(M0) - 10.7;
    e = floor(log10(max(abs(M))));
    ev(k).fac = 10^(e - 20);                     % idagrn: M = mt * fac * 1e20
    ev(k).mt = M / (1e20 * ev(k).fac);
end
end
