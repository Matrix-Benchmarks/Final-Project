function [X_hat,Out] = MS_FGD(U0,f_grad,eta,~,T,~,max_time)
%n=length(A);
% Initialization1
% if init==1
% [U0,V0] = initialization_ms(y,n,r,d1,d2,Omega);
% end
% if init==2
% % Initialization2
% [U0,V0] = initialization_gdms(y,A,n,r,tau,initialeta,rd);
% end
% if init==3
% % Initialization3
% [U0,V0] = initialization_rdms(y,A,n,r,0,0.015);
% end
%loss = zeros(T,1);
%dist = zeros(T,1);
%time = zeros(1,T);
%X_hat = cell(1,T);
% Projected GD
U=U0;
%V=V0;
tic;
%dist(1) = norm(U*U'- X_star, 'fro');
t=1;
while toc < max_time
    % Calculate the gradient
    %
    f_gradX = f_grad(U * U');
    nabla_U = 2*f_gradX * U; 
    
    %update
    U = U - eta * nabla_U;
    
    
    %dist(t) = norm(U*U'- X_star, 'fro');
    time(t) = toc;
    X_hat{t} = {U,U};
    
    t = t+1;
end
toc;
X_hat = X_hat(1:t-1);
Out.time = time(1:t-1);
end
