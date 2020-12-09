function [M_star] = data_generation_ms(d1,d2,r)

% Generate M* 
U_star = randn(d1,r) / sqrt(d1);   % randn() is a little more efficent (to normrnd() )
M_star = U_star*U_star';
end