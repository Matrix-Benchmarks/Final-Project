% =========================================================================
%             Matrix Algebraic Pursuit algorithms - Demo
% =========================================================================
% Matrix Algebraic Pursuit (ALPS) algorithms are accelerated hard thresholding 
% methods on low rank recovery of linear inverse systems. In this example, we
% solve the problem where:
%                         y = A(X*) + e
%             A   : m x n -> p x 1 mapping 
%             X*  : m x n rank-k data matrix
%             e   : p x 1 additive noise vector 
% Matrix ALPS solve the following minimization problem 
%          minimize ||y - A(X)||_F^2    subject to rank(X) <= k. 
% 
% Detailed discussion on the algorithm can be found in 
% [1] "Matrix Recipes for Hard Thresholding Methods", written
% by Anastasios Kyrillidis and Volkan Cevher, Technical Report, 2011.
% =========================================================================
% 14/01/2012, by Anastasios Kyrillidis. anastasios.kyrillidis@epfl.ch, EPFL.
% =========================================================================

close all
clear
clc

global global_Mat;

addpath ALPS/;
addpath operators/;
addpath utility/;
addpath PROPACK/;

n = 2000;                       % Number of columns
m = 800;                        % Number of rows
k = 20;                         % Rank of X*
p = floor(0.15*n*m);            % Number of observations
sigma = 0*10^-4;                % Noise variance

global_Mat = zeros(m, n);

% ALPS parameters
params.ALPSiters = 500;         % Maximum number of iterations
params.tol = 5e-6;              % Convergence tolerance
params.xpath = 1;               % Keep history log
params.svdMode = 'propack';     % SVD method - default mode for Matrix ALPS II - other options: 'svds', 'svd'
params.cg_tol = 1e-10;          % Conjugate gradients tolerance
params.cg_maxiter = 500;        % Maximum number of conjugate gradient iterations
params.svdApprox = 0;           % Set to 1 for column subset selection - really slow...
params.power = 2;               % Number of iterations in subspace range finder

maxits = 0;
MC = 1;                         % Monte Carlo iterations
for i = 1:MC
    L = randn(m, k); R = randn(k, n);   % Left and right component such that X* = L*R;
    X = L*R;                            
    X = X./norm(X, 'fro');              % Normalize such that ||X*||_F = 1
    Xt = X';                            % Transpose of X' -> used in QR method

    %% Generate linear mapping operator - taken from SparCS implementation

    % Options: prob = 'noiselets' -----> Noiselets as linear operator - IMPORTANT: use power of 2 for matrix size
    %          prob = 'robustMC' ------> Robust Matrix Completion setting

    prob = 'robustMC';
    switch prob
        case 'noiselets'
            % Noiselets as linear operator A                    
            [ix iy p] = randmask(m, n, p);
            linearInd = sub2ind([m n], double(ix), double(iy)); 
            tlinearInd = sub2ind([n m], double(iy), double(ix)); 

            linearInd = linearInd(:);
            tlinearInd = tlinearInd(:);

            idx = linearInd;
            tidx = tlinearInd;
            idx2 = randperm(n*m);
            
            A = @(z) Anoiselet(z,idx, idx2);
            At = @(z) Atnoiselet(z,idx, idx2, size(X));    

            tA = @(z) Anoiselet(z, tidx, idx2);
            tAt = @(z) Atnoiselet(z,tidx, idx2, size(Xt));            
            
        case 'robustMC'
            % Robust matrix completion
            [ix iy p] = randmask(m, n, p);
            linearInd = sub2ind([m n], double(ix), double(iy)); 
            tlinearInd = sub2ind([n m], double(iy), double(ix)); 

            linearInd = linearInd(:);
            tlinearInd = tlinearInd(:);

            idx = linearInd;
            tidx = tlinearInd;

            A = @(z) z(idx);
            At = @(z) full(sparse(double(ix), double(iy), z, m, n));

            tA = @(z) z(tidx);
            tAt = @(z) full(sparse(double(iy), double(ix), z, n, m));        
    end 

    %% Generate noise and observations
    e = randn(m*n,1);                   
    e_sparse = zeros(m*n,1);
    e_sparse(idx) = e(idx);
    e = e_sparse;
    e = sigma*e/norm(e);
    e_sparse = e;

    if (sigma == 0)
        y = A(X);
        yt = tA(Xt);
    else
        e_sparse(e_sparse == 0) = [];
        y = A(X) + e_sparse;
        yt = tA(Xt) + e_sparse;
    end;

    %% Matrix ALPS I
    disp('=====================================================================');
    t0 = clock;
    [X_hat, numiter, X_path] = matrixALPSI(y, A, At, m, n, k, params, X);
    toc_time = etime(clock, t0);
    if (numiter ~= -1)
        str = sprintf('MatrixALPS I terminated - Error norm: %f ', norm(X-X_hat, 'fro')/norm(X, 'fro'));
        disp(str);
        str = sprintf('Number of iterations: %f ', numiter);
        disp(str);        
    end;

    if (numiter > maxits)
        maxits = numiter;
    end;

    num_iter(1, i) = numiter;
    error(1,i,1:numiter) = X_path;
    time(1, i) = toc_time;
    rec_error(1,i) = norm(X-X_hat, 'fro')/norm(X, 'fro');

    global_Mat = zeros(m, n);
    %% Matrix ALPS II
    params.tau = 0;
    disp('=====================================================================');
    t0 = clock;
    [X_hat2, numiter2, X_path2] = matrixALPSII(y, A, At, m, n, k, params, X);
    toc_time = etime(clock, t0);
    if (numiter2 ~= -1)
        str = sprintf('MatrixALPS II terminated - Error norm: %f ', norm(X-X_hat2, 'fro')/norm(X, 'fro'));
        disp(str);
        str = sprintf('Number of iterations: %f ', numiter2);
        disp(str);        
    end;

    if (numiter2 > maxits)
        maxits = numiter2;
    end;

    num_iter(2, i) = numiter2;
    error(2,i,1:numiter2) = X_path2;
    time(2, i) = toc_time;
    rec_error(2,i) = norm(X-X_hat2, 'fro')/norm(X, 'fro');

    global_Mat = zeros(m, n);
    %% Matrix ALPS II with randomized subspace finder
    params.tau = 0;
    disp('=====================================================================');
    t0 = clock;
    [X_hat3, numiter3, X_path3] = matrixALPSII_QR(yt, tA, tAt, n, m, k, params, Xt);
    toc_time = etime(clock, t0);
    if (numiter3 ~= -1)
        str = sprintf('MatrixALPS II QR terminated - Error norm: %f ', norm(Xt-X_hat3, 'fro')/norm(Xt, 'fro'));
        disp(str);
        str = sprintf('Number of iterations: %f ', numiter3);
        disp(str);
    end;

    if (numiter3 > maxits)
        maxits = numiter3;
    end;

    num_iter(3, i) = numiter3;
    error(3,i,1:numiter3) = X_path3;
    time(3, i) = toc_time;
    rec_error(3,i) = norm(Xt-X_hat3, 'fro')/norm(Xt, 'fro');
    
    global_Mat = zeros(m, n);
    %% Matrix ALPS II with randomized subspace finder and orthogonalization
    params.tau = 0;
    disp('=====================================================================');
    t0 = clock;
    [X_hat4, numiter4, X_path4] = matrixALPSII_QR_ortho(yt, tA, tAt, n, m, k, params, Xt);
    toc_time = etime(clock, t0);
    if (numiter4 ~= -1)
        str = sprintf('MatrixALPS II QR and ortho terminated - Error norm: %f ', norm(Xt-X_hat4, 'fro')/norm(Xt, 'fro'));
        disp(str);
        str = sprintf('Number of iterations: %f ', numiter4);
        disp(str);
    end;

    if (numiter4 > maxits)
        maxits = numiter4;
    end;

    num_iter(4, i) = numiter4;
    error(4,i,1:numiter4) = X_path4;
    time(4, i) = toc_time;
    rec_error(4,i) = norm(Xt-X_hat4, 'fro')/norm(Xt, 'fro');
end;

set(0, 'DefaultAxesFontSize', 16);
figure;
semilogy(squeeze(median(error(1,:,:),2)),'k','LineWidth',4); hold on
semilogy(squeeze(median(error(2,:,:),2)),'g','LineWidth',4);  hold on
semilogy(squeeze(median(error(3,:,:),2)),'r--','LineWidth',4); hold on
semilogy(squeeze(median(error(4,:,:),2)),'b--','LineWidth',4); hold on
grid on;
xlabel('\# of iterations','FontSize', 18, 'Interpreter', 'latex')
ylabel(strcat(['$\|\mathbf{X}(i) $','- ','$ \mathbf{X}^\ast\|_F$']),'Interpreter', 'latex','FontSize', 18)
axis([1 maxits + 20 10^-5 1]);
semilogy((norm(e))*ones(1,maxits + 20),'k-.','LineWidth',2);

legend(strcat('Matrix ALPS I [', num2str(median(time(1,:))),']'), ...        
       strcat('Matrix ALPS II [', num2str(median(time(2,:))),']'), ...;
       strcat('Matrix ALPS II with QR [', num2str(median(time(3,:))),']'), ...;
       strcat('Matrix ALPS II with QR and ortho [', num2str(median(time(4,:))),']'));   
h  = legend;
set(h, 'interpreter', 'latex','FontSize', 18);
shg;

    
    