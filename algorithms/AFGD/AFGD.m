function [Mout, Out] = AFGD(m,n,r,Omega,data,start,max_time)
%
% Implementation of matrix completion algorithm, Accelerated Factored
% Gradient Descent
% Code by Richard Morse
% =========================================================================
% Adapted from algorithm in "Accelerated Factored Gradient Descent for
% Low-Rank Matrix Factorization" by Dngruo Zhou, Yuan Cao, and Quanquan Gu
% =========================================================================
% Set parameters
maxiter = 1e7; % Very large number just to make sure we never reach it
saveiterates = 1;
itres = zeros(maxiter,1);
conv_rate = 0;
p = length(data);
verbosity = 1;

% creat a sparse matrix
[I, J] = ind2sub([m,n],Omega);
diff_on_omega_matrix = sparse(I,J,data,m,n,p); %data vector in sparse matrix format

if ~isempty(start)
  U0 = start.U;
  if verbosity
    fprintf('initial point provided!\n')
  end
else
  U0 = randn(m,n);
  if verbosity
    fprintf('no initial point!\n')
  end
end
V = U0;
X = U0;


clear opts;
clear start;

[A0,D0,B0] = svds(X,r);
%for ACCPROJ algorithm
D0inv = inv(D0);
inner_eta = D0(r,r)^2; 
inner_beta = (norm(D0,2) - D0(r,r)) / (norm(D0,2) + D0(r,r));

eta = 0.005;
gamma = 0.05;%D0(r,r)^6 / norm(U0)^4;
alpha = sqrt(eta*gamma);
inner_maxiter = 10;

%initial Y
Y = (alpha/(alpha+1))*V + (1/(alpha+1))*X;
%gradient
diff_on_omega = data' - partXY(Y',Y',I,J,p); %derivative of obj at initial point
updateSval(diff_on_omega_matrix,diff_on_omega,p); %putting derivative into sparse matrix format
res = norm(diff_on_omega);

iter = 1;
itres(iter) = res;
time = zeros(1, maxiter);
tic

% iteration
while  toc < max_time
    disp(iter)
    %updates
    Vp = (1-alpha)*V + alpha*Y - (alpha/gamma)*diff_on_omega_matrix*Y;
    V = accproj(Vp,A0,B0,D0inv,r,inner_eta,inner_beta,inner_maxiter);
    
    [P,~,Q] = svd((Y - eta*diff_on_omega_matrix*Y)'*U0);
    X = V*(P*Q'); %projecting X onto U0 PSD space
    
    %doing this last rather tahn first to make iterations line up with other algos
    Y = (alpha/(alpha+1))*V + (1/(alpha+1))*X;
    %gradient
    diff_on_omega = data' - partXY(Y',Y',I,J,p); %derivative of obj at initial point
    updateSval(diff_on_omega_matrix,diff_on_omega,p); %putting derivative into sparse matrix format
    %residual - norm of gradient which is also error b/w Y and U* cuz grad_f(U) = -(UU' - data)
    %calculating on Y
    res = norm(diff_on_omega);
    
    if saveiterates
        Mout{iter}={X,X};
    end

    itres(iter) = res;
    conv_rate = (itres(iter)/itres(max(1,iter-15)))^(1/min(15,iter-1));
    time(iter) = toc;
    iter = iter + 1;
end
iter=iter-1;
if saveiterates
    Mout = Mout(1:iter);
else
    Mout{iter}={X,X'};
end

Out.itrelres = itres(1:iter)/norm(data);
Out.iter = iter;
Out.reschg = abs(1-conv_rate);
Out.time = time(1:iter);
end