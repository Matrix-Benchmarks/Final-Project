function [loss] = loss_ms(U,V,A,y)
X=U*V';
n=length(A);
loss=0;
for i=1:n
    loss=loss+(sum(sum(A(:,:,i).*(X)))-y(i))^2;
end
loss=loss/(2*n)+1/8*norm(U'*U-V'*V,'fro')^2;
