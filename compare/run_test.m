function [Xr,outs] = run_test(Phi,y,r,alg_names, max_time, output_file, goal_matrix, filter)
% This function runs different algorithms (indicated by alg_name)
% for a given matrix completion problem with entry mask Phi and provided 
% data y up to a max time.
% =========================================================================
% Parameters
% ----------
% Phi:  (d1 x d2) sparse matrix containing m non zero entries.  The indices
% of the non zero entries correspond to the data vector (y).
% y:    (m x 1) vector of samples. Provided matrix entry values,
%       ordered in according to linear indices of non-zero entries of Phi.
% r:    Target rank to be used.
% alg_name:     (1 x nr_algos) cell of character strings. Indicating which 
%               algorithm to use.
% Returns
% ----------
% Xr:   (1 x nr_algos_total) cell. The i-th cell contains algorithmic
%       iterates of i-th algorithm.
% outs: (1 x nr_algos_total) cell. The i-th cell contains additional
%       information about the progress of the i-th algorithm.
% Author: Josh Engels, adopted from test code by Christian Kuemmerle 2020.

[d1,d2]=size(Phi);
Omega = find(Phi);
[I, J] = ind2sub([d1 d2],Omega);
nr_algos = length(alg_names);

for alg_num = 1:nr_algos
    current_alg = alg_names{alg_num};
    
    if (any(strcmp(["ScaledASD", "ASD", "NIHT_Matrix", "CGIHT_Matrix"], current_alg)))
        alg_func = str2func(current_alg + "_outp");
        start = make_start_x_IHT_ASD(current_alg,d1,d2,r,Omega,y);
        [iterations,times] = alg_func(d1,d2,r,Omega,y,start,max_time);
    end
    
    if(current_alg == "R3MC") 
        [rowind,colind]=find(Phi);
        input = struct;
        input.d1     = d1;
        input.d2     = d2;
        input.r      = r;
        input.Phi    = Phi;
        input.y      = y;
        input.data_ls.rows = rowind;
        input.data_ls.cols = colind;
        input.data_ls.entries = input.y;
        input.data_ls.nentries = length(input.y);
        input = initialize_R3MC(input,false);
        [~, times] = R3MC_adp(input, max_time);
        iterations = times.X;
    end
    
    if strcmp(current_alg,'LMaFit')
        [~,~,Out] = lmafit_mc_adp(d1,d2,r,Omega,y, max_time);
        times=Out;
        times.N=Out.iter;
        iterations = cell(1,times.N);
        for kk=1:times.N
            iterations{kk}={Out.Xhist{kk},Out.Yhist{kk}'};
        end
    end
    
    if strcmp(current_alg,'ScaledGD')
        input.d1 = d1;
        input.d2 = d2;
        [input.rowind,input.colind]=find(Phi);
        [input.Omega] = find(Phi);
        input.y = y;
        [iterations,times] = ScaledGD(input, r, max_time);
    end
    
    if strcmp(current_alg, 'MatrixIRLS')        
        opts = getDefaultOpts_IRLS;
        input.d1     = d1;
        input.d2     = d2;
        input.r      = r;
        input.Phi    = Phi;
        input.y      = y;
        [~,times] = MatrixIRLS(input,0, opts, max_time);
        iterations = cell(1,times.N);
        for kk=1:times.N 
            iterations{kk}=times.X{kk};
        end
    end
    
    if (any(strcmp(["FGD", "AGD", "AFGD"], current_alg)))
        alg_func = str2func("MS_" + current_alg);
        [start,~] = initialization_ms(y,r,d1,d2,Omega);
        num_observations = length(y);
        
        At_id = @(z) full(sparse(I,J,z,d1,d2,num_observations));
        f_grad = @(x) (-1)*At_id(y - x(Omega));
        
        [~,Sigma1,~] = svds(goal_matrix,r); Sigma1 = diag(Sigma1); sigma1 = Sigma1(1);
        eta = 1;%0.2/sigma1;
        gamma = 1;%(0.1^2)/eta;
        
        T=10000;
        TS = 70;
        
        [iterations,times,dist] = alg_func(start,goal_matrix,f_grad,eta,gamma,T,TS,max_time);
        disp(dist(end))
    end
    
    if strcmp(current_alg, 'MatrixALPS')
        params.maxtime = max_time;       % Max time for algorithm to run
        params.ALPSiters = 500;         % Maximum number of iterations
        params.tol = 5e-12;             % Convergence tolerance (not used)
        params.xpath = 1;               % Keep history log
        params.svdMode = 'propack';     % SVD method - default mode for Matrix ALPS II - other options: 'svds', 'svd'
        params.cg_tol = 1e-10;          % Conjugate gradients tolerance
        params.cg_maxiter = 500;        % Maximum number of conjugate gradient iterations
        params.svdApprox = 0;           % Set to 1 for column subset selection - really slow...
        params.power = 2;               % Number of iterations in subspace range finder
        
        A = @(z) z(Omega);
        At = @(z) full(sparse(I, J, z, d1, d2));
        
        %Only ALPS I & II integrated so far (neither QR variation has been
        %implemented)
        [~,times,iterations] = matrixALPSII(y, A, At, d1, d2, r, params, goal_matrix);
    end
    
    print_result_to_file(current_alg, times, iterations, output_file, goal_matrix, filter);
    
end