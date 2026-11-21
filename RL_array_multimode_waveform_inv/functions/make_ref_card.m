function card = make_ref_card(mdir, refid, premid, trueid, zmoho)
% card = make_ref_card(mdir, refid, premid, trueid, zmoho)
%
% Reference ("input") model for the recovery test: PREM below the Moho of
% the true model, and the true model's water, sediments and crust above it
% (the shallow structure at NoMelt is known from bathymetry and active-
% source data; a 24-km PREM crust under 3 km of water would bias any mantle
% inversion at an oceanic site). Between the Moho and PREM's own Moho
% (24.4 km) the PREM uppermost-mantle values are extended upward.
%
% Writes <mdir>/CARDS/<refid>.card and copies <premid>.qmod to <refid>.qmod.
%
% mdir   - path to run_MINEOS
% refid  - name of the reference card to write (e.g. 'prem_nomelt')
% premid - PREM card (e.g. 'prem_35')
% trueid - card whose water/crust is used (e.g. 'pa5_5km')
% zmoho  - Moho depth below sea level (km) in the true card (PA5: 12)
%
% jbrussell 11/2026

cp = read_model_card(fullfile(mdir, 'CARDS', [premid, '.card']));
ct = read_model_card(fullfile(mdir, 'CARDS', [trueid, '.card']));
rmoho = (6371 - zmoho)*1000;
f = {'rad','rho','vpv','vsv','qkap','qmu','vph','vsh','eta'};

% PREM mantle: nodes below PREM's Moho (deepest node of each PREM node
% pair at the Moho radius is the mantle side)
rpm = max(cp.rad(cp.vsv > 4000 & cp.rad < 6360000));        % PREM Moho (mantle side)
ip = find(cp.rad <= rpm);
ip = ip(1:find(cp.rad(ip) == rpm, 1, 'first'));            % mantle-side node at the PREM Moho
% one node at the true Moho with the PREM uppermost-mantle values
it = find(ct.rad >= rmoho);
it = it(find(ct.rad(it) == rmoho, 1, 'last'):end);           % crust side at the Moho and above
for k = 1:numel(f)
    v = cp.(f{k})(ip);
    top = cp.(f{k})(ip(end)); if strcmp(f{k}, 'rad'), top = rmoho; end
    card.(f{k}) = [v; top; ct.(f{k})(it)];
end
card.z = 6371 - card.rad/1000;
write_MINEOS_mod(card, fullfile(mdir, 'CARDS', [refid, '.card']));
copyfile(fullfile(mdir, 'CARDS', [premid, '.qmod']), fullfile(mdir, 'CARDS', [refid, '.qmod']));
end
