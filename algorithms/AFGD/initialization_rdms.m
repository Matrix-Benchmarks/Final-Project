function [U0,V0] = initialization_rdms(y,A,n,r,mu,nu)
[d1,d2] = size(A(:,:,1));
%[d1,d2] = size(X);
X=zeros(d1,d2);
for i=1:n
X=X+y(i)*A(:,:,i);
end
[U,Sigma,V] = svds(X/n,r);
Sigma = sqrt(Sigma);
U0=U*Sigma+normrnd(mu,nu,d1,r);
V0=V*Sigma+normrnd(mu,nu,d2,r);
end