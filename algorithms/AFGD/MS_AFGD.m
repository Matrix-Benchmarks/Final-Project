function [X_hat,Out,dist] = MS_AFGD(U0,X_star,f_grad,eta,gamma,T,TS,max_time)
%Inputs
%U0         : initial guess
%X_star     : ground truh matrix if computing error locally
%f,f_grad   : objective function and its gradient
% eta,gamma : step_size, acceleration
%T          : max # of iterations
%TS         : max # of iterations for internal proj algorithm
%max_time   : max time alloted to run algorithm for
%currently only using initialization scheme 1, initialization_ms

%Outputs
%X_hat : solution iteratives, saved as a cell array of the factor {U,U}
%Out   : solution information, only saving time information as of now
%dist  : distance function if computed locally, used within AFGD package,
%        but not within framework
%======================================================================
%Oringally by Dr. Dongruo Zhou, updated by Ricky More, 10 Decmeber 2020
%======================================================================
% n=length(A);
% % Initialization1
% if init==1
% [U0,V0] = initialization_ms(y,A,n,r);
% end
% if init==2
% % Initialization2
% [U0,V0] = initialization_gdms(y,A,n,r,tau,initialeta,rd);
% end
% if init==3
% % Initialization3
% [U0,V0] = initialization_rdms(y,A,n,r,0,0.015);
% end

%dist = zeros(1,T);
%time = zeros(1,T);
%X_hat = cell(1,T);
% Projected GD
U=U0;
[S1,D1,S2] = svd(U0, 'econ');
V=U;
W= U;
alpha = sqrt(eta*gamma);
tic;
dist(1) = norm(U*U'- X_star, 'fro');
t=1;
while toc < max_time
    % Calculate the gradient
    %
    U =  alpha/(alpha+1)*V + 1/(alpha+1)*W;
    f_gradX = f_grad(U * U');
    nabla_U = 2*f_gradX * U; 
    
    %update
    V = (1-alpha)*V + alpha*U - alpha/gamma*nabla_U;
    V = proj(V, S1,D1,S2,TS);
    W = rotation(U - eta * nabla_U, U0);

    
    dist(t) = norm(W*W'- X_star, 'fro');
    time(t) = toc;
    X_hat{t} = {W,W};
    t = t+1;
end
toc;
X_hat = X_hat(1:t-1);
Out.time = time(1:t-1);
end
