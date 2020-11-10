%  TODO: Understand parameters

function  opts_custom = default_params(cond_nr, m, n)

% %%% All parameters should use a time limt
opts_custom.time_limit_ms = 1e6;



% %%% Optional parameters addressing options of 'second-order algorithms'.
opts_custom.N0 = 400;           % Max. number of (outer) iterations for 
                                % 'second-order algorithms', which include 
                                % MatrixIRLS, R2RILS and RTRMC.
opts_custom.tol_CG_fac=1e-5*cond_nr^(-1);    % tolerance for stopping criterion of inner iterations
opts_custom.N0_inner = 500;     % Max. number of (inner) iterations for 'second-order algorithms'

%%% Optional parameters addressing only options for IRLS
opts_custom.p = [0]; % (non-)convexity parameters p for IRLS
% p = 0: sum of log objective
% 0 < p < 1: Schatten-p quasi-norm.
% p = 1: Nuclear norm.
opts_custom.R = min(m,n);             %min(d1,d2);%floor(10*r);50;%
opts_custom.adaptive_cg = 0;
opts_custom.mode_linsolve = 'tangspace';    % 'rangespace'
opts_custom.type_mean = {'optimal'};        % 'harmonic','arithmetic','geometric','min'
opts_custom.objective = 'objective_thesis'; % 'pluseps','pluseps_squared'
opts_custom.mode_eps = 'oracle_model_order';% 'iter_diff','auto_decay'
opts_custom.epsmin = 1e-16;
opts_custom.tracking = 0;
opts_custom.lambda = 0;
opts_custom.increase_antisymmetricweights=0;
opts_custom.saveiterates = 1;
opts_custom.verbose = 1;
