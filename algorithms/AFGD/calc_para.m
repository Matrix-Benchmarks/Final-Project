function [mu,sigma1,cond_num] = calc_para(M_star,r)

[d1,d2] = size(M_star);
[A,Sigma,B] = svds(M_star,r);  % Singular Value decomposition
sigma = diag(Sigma);

% Get the parameter of M_star

a1 = max(sum(A.^2,2));
a2 = max(sum(B.^2,2));
% temp = A*B';
% a3 = max(temp(:));

% Calculate the minimum value of incoherence mu
mu = max(a1*d1/r,a2*d2/r);
% mu = max([a1*d/r,a2*d/r,a3^2*d^2/r]);

% Calculate matrix 2 norm
sigma1 = sigma(1);


% Calculate condition number
cond_num = sigma(1)/sigma(r);

end