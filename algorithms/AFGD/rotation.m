function [yy] = rotation (X,Y)
%%%Y is the base

[A,~,B] = svd(X'*Y);
yy = X*A*B';




end