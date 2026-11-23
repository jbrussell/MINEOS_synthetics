function parameters = aswms_setup(opt, component)
% parameters = aswms_setup(opt, component)
% Build the ASWMS-Q parameters for one component (aswms_parameters.m) and
% install them for the ASWMS-Q functions in functions/aswms, which read them
% through setup_parameters.m.
global ASWMS_PARAMETERS_
parameters = aswms_parameters(opt, component);
ASWMS_PARAMETERS_ = parameters;
end
