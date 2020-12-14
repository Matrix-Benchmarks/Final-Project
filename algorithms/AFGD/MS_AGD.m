function [X_hat,Out,dist] = MS_AGD(U0,X_star,f_grad,eta,gamma,T,~,max_time)
%======================================================================
%Oringally by Dr. Dongruo Zhou, updated by Ricky More, 10 Decmeber 2020
%======================================================================
% n=length(A);
% %Initialization1
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
%dist = zeros(T,1);
%X_hat = cell(T,1);
%time = zeros(T,1);
% Projected GD
U=U0;
V=U;
W= U ;
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
    W = U - eta * nabla_U;
    V = (1-alpha)*V + alpha*U - alpha/gamma*nabla_U;
    
    
    dist(t) = norm(W*W'- X_star, 'fro');
    time(t) = toc;
    X_hat{t} = {W,W};
    
    t = t+1;
end
toc;
X_hat = X_hat(1:t-1);
Out.time = time(1:t-1);
end
