function [U0,V0] = initialization_gdms(y,A,n,r,tau,initialeta,rd)
[d1,d2] = size(A(:,:,1));
%[d1,d2] = size(X);
if rd==0
X=zeros(d1,d2);
for i=1:n
X=X+y(i)*A(:,:,i);
end
[U,Sigma,V] = svds(X/n,r);
Sigma = sqrt(Sigma);
U0=U*Sigma;
V0=V*Sigma;
X0=U0*V0';
end
if rd==1
X=normrnd(0,1,d1,d2);
[U,Sigma,V] = svds(X/n,r);
Sigma = sqrt(Sigma);
U0=U*Sigma;
V0=V*Sigma;
X0=U0*V0';
end
mat=zeros(d1,d2,n);
%gd initialization
%step size
eta=initialeta;
for t = 1:tau 
    % Calculate the gradient
    for i=1:n
    mat(:,:,i)=(sum(sum(A(:,:,i).*(X0)))-y(i))*A(:,:,i);
    end
    %
    nablaX = 1/n*sum(mat(:,:,:),3);
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