function start = make_start_x_AFGD(func_name,m,n,r,Omega,data)
%pulled from Theorem 11 in http://proceedings.mlr.press/v49/bhojanapalli16.pdf
%as recomended by appendix of
%http://proceedings.mlr.press/v108/zhou20b/zhou20b.pdf (AFGD paper)

% creat a sparse matrix
[I, J] = ind2sub([m,n],Omega);
p = length(Omega);
M_omega = sparse(I,J,data,m,n,p);

%gradient at 0 is -M_omega
%projecting gradient at 0 onto PSD cone
[V,D] = eigs(M_omega,m); grad0PSD = V * max(D,0) * V';

X0 = 1*grad0PSD; %term in front ends up being 1= 1/1 since bottom simplifies to norm
% of e1*e1' in matrix completion case, which is obs just 1

[V,D] = eigs(X0,r); %X0 = VDV' w/ V rank-r
U0 = V*sqrt(D); %U0 = Vsqrt(D) - lol no idea if this makes any sense
start.U = U0;
%start.U = zeros(m,r); %might be better

end