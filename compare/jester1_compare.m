% Similar to protocal described here
% http://proceedings.mlr.press/v70/wang17n/wang17n.pdf

% Read in Jester1 data
file_name = "jester1.csv";
jester_data = readmatrix("datasets/" + file_name);

% Get dimensions of data
num_rows = size(jester_data, 1);
num_cols = size(jester_data, 2);

% Guess for the rank
% TODO: Try each algo with different guesses? Sensitivty to the "wrong"
% rank, more realistic to the real world?
r = 5;

% Say number of measurement = 1 / 20th of the matrix, can vary for different tests
m = cast(num_cols * num_rows / 2, 'uint32');

% Get the mask
% Phi is a sparse matrix representing a mask of an arbitrary d1*d2 matrix,
% where an entry is 1 if it should be included and 0 otherwise.
% Phi contains m 1s.
% resample means that the sampling procedure will keep trying (up to a
% max of max_nr_resample times) to get r entries per row and per column.
% Omega is a linear array of the indices of the 1s in Phi
max_nr_resample = 100;
[Phi,Omega] = sample_phi_MatrixCompletion(num_rows,num_cols,m,'resample',r,max_nr_resample);

% Linear array of observations
observations = jester_data(Omega);


%% Temporary
%% Choose algorithms for matrix completion
% alg_names={'MatrixIRLS','R3MC','R3MC-rankupd','R2RILS','RTRMC','LRGeomCG',...
%     'LMaFit','ScaledASD','ScaledGD','NIHT'};
alg_names = {'ScaledASD'};
% 
% %%% Set optional algorithmic parameters
opts_custom.tol = 1e-10;        % tolerance for stopping criterion
opts_custom.N0 = 400;           % Max. number of (outer) iterations for 
                                % 'second-order algorithms', which include 
                                % MatrixIRLS, R2RILS and RTRMC.
opts_custom.N0_firstorder = 1000; % Max. number of iterations for 'first-order algorithms'.

% %%% Optional parameters addressing options of 'second-order algorithms'.
opts_custom.tol_CG_fac=1e-5*cond_nr^(-1);    % tolerance for stopping criterion of inner iterations
opts_custom.N0_inner = 500;     % Max. number of (inner) iterations for 'second-order algorithms'


%%% Optional parameters addressing only options for IRLS
opts_custom.p = [0]; % (non-)convexity parameters p for IRLS
% p = 0: sum of log objective
% 0 < p < 1: Schatten-p quasi-norm.
% p = 1: Nuclear norm.
opts_custom.R = min(d1,d2);             %min(d1,d2);%floor(10*r);50;%
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

%% Run algorithms for matrix completion
[Xr,outs,alg_names] = run_MC_algos(Phi,observations,r,alg_names,opts_custom);
