% D2_STEP2_PREP_MINEOS_IDAGRN  Step 2, part 1 (driver): MINEOS and idagrn6.
%
% With the step-1 model (CARDS/<opt.s1_final>.card, from D1) as the
% reference model (or, opt.ref_card_id set, an existing card: e.g. the true
% model opt.trueid for an oracle test), runs once and caches in opt.prepfile:
%   - MINEOS mode tables of the true model (PA5, branches
%     0..opt.nbranch_data-1) and of the reference (branches 0..opt.nm-1)
%   - reference phase/group velocities and Vsv, Vsh, Vpv, Vph, rho kernels
%     of branches 0..opt.nm-1 at the dispersion nodes
%   - idagrn6 seismograms at the 16 NoMelt OBSs for the 15 events:
%       data      : opt.data_source = 'aswms' (default): NOT computed here;
%                   the data are the ASWMS-Q event files of A1-A3 (for the
%                   synthetic test from the SAC files of A0). 'internal':
%                   PA5, all opt.nbranch_data branches (T and Z)
%       Uref      : reference model, one branch at a time (idagrn6_mask),
%                   also written in the ASWMS-Q event-file format to
%                   <opt.aswms_dir>branchmat/<refid>/<event id>_<component>.mat
%                   (write_branchmat.m; stadata(j).data = [Nt x nm])
% (functions/prep_mineos_idagrn.m). Takes a few minutes.
% Then record sections of the data and of the step-1 model (T and Z,
% inversion band, with the automatic windows of A3 when they have been
% picked) in figs_step2/ (functions/plot_record_sections.m).
%
% Needs the MINEOS and idagrn6 binaries (../FORTRAN/bin) and SAC (sacio).
% Previous: D1_step1_rayleigh_isotropic_inv.m and (data) A1-A3 (synthetic
% test: A0 first). Next: D3_step2_waveform_inv.m

parameters_RL;
o = opt;
if isempty(opt.ref_card_id)
    S1 = load(opt.s1_file, 'card_final');
    o.ref_card = S1.card_final;             % the step-1 isotropic model
    o.refid = opt.s1_final;
else                                        % an existing card (e.g. the true model: oracle test)
    o.ref_card = opt.ref_card_id; o.refid = opt.ref_card_id;
    fprintf('reference model: CARDS/%s.card (opt.ref_card_id; not the step-1 model)\n', opt.ref_card_id);
end
P = prep_mineos_idagrn(o);
fprintf('cached in %s: data branches 0-%d, modeled branches 0-%d\n', opt.prepfile, opt.nbranch_data-1, opt.nm-1);
Pp = P;                                     % record sections of the data: from the event files
win_ = {{}, {}};                            % automatic windows of A3 (event.winpara), if picked yet
ty_ = 'TS';
for wv = 1:2
    for e = 1:numel(P.ev)
        Dd = [];
        if strcmp(opt.data_source, 'aswms') || exist(fullfile(opt.aswms_dir, 'eventmat'), 'dir') == 7
            Dd = read_event_data(opt, P.ev(e), P.sta.name, opt.aswms_components{wv}, P.t);
        end
        if strcmp(opt.data_source, 'aswms')
            if isempty(Dd), u = zeros(numel(P.t), numel(P.sta.name)); else, u = Dd.U; u(~isfinite(u)) = 0; end
            Pp.(['u_true_' ty_(wv)]){e} = u;
        end
        if ~isempty(Dd) && ~isempty(Dd.winpara)
            win_{wv}{e} = struct('winpara', Dd.winpara, 'dist', Dd.dist, 'toff', Dd.toff);
        else
            win_{wv}{e} = [];
        end
    end
end
plot_record_sections(Pp, 'figs_step2', optget(opt, 'fc', []), win_);
fprintf('record sections in figs_step2/\n');
