function [U0,V0] = initialization_rdms(y,r,d1,d2,Omega,mu,nu)
%======================================================================
%Oringally by Dr. Dongruo Zhou, updated by Ricky More, 10 Decmeber 2020
%======================================================================
[I, J] = ind2sub([d1 d2],Omega);
X = sparse(I,J,y,d1,d2,length(y));

[U,Sigma,V] = svds(X,r);
Sigma = sqrt(Sigma);
U0=U*Sigma+normrnd(mu,nu,d1,r);
V0=V*Sigma+normrnd(mu,nu,d2,r);
end