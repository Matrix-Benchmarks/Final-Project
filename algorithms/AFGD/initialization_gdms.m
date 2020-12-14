function [U0,V0] = initialization_gdms(y,r,d1,d2,Omega,tau,initialeta,rd)
%======================================================================
%Oringally by Dr. Dongruo Zhou, updated by Ricky More, 10 Decmeber 2020
%======================================================================
% if rd==0
% [I, J] = ind2sub([d1 d2],Omega);
% X = sparse(I,J,y,d1,d2,length(y));
% [U,Sigma,V] = svds(X/n,r);
% Sigma = sqrt(Sigma);
% U0=U*Sigma;
% V0=V*Sigma;
% X0=U0*V0';
% end

if rd==1
X=normrnd(0,1,d1,d2);
[U,Sigma,V] = svds(X,r);
Sigma = sqrt(Sigma);
U0=U*Sigma;
V0=V*Sigma;
X0=U0*V0';
end
%gd initialization
%step size
eta=initialeta;
[I, J] = ind2sub([d1 d2],Omega);
for t = 1:tau 
    % Calculate the gradient
    mat = sparse(I,J,X0(Omega) - y,d1,d2,length(y));
    nablaX = mat;
    %update
    X0=X0-eta*nablaX;
    %rank-r projection
   [U,Sigma,V] = svds(X0,r);
   Sigma = sqrt(Sigma);
   U0=U*Sigma;
   V0=V*Sigma;
   X0=U0*V0';
end


end