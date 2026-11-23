% SETUP_PARAMETERS  (this package's replacement for setup_parameters.m of
% ASWMS-Q) The ASWMS-Q functions in this folder call setup_parameters to get
% the "parameters" structure. Here it returns the structure installed by
% aswms_setup.m (from opt, parameters_RL.m, for the current component), or
% the package defaults if none is installed.
global ASWMS_PARAMETERS_
if isempty(ASWMS_PARAMETERS_)
    parameters = aswms_parameters(struct(), 'BHZ');
else
    parameters = ASWMS_PARAMETERS_;
end
