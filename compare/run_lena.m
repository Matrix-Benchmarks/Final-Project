%% Intro
% Adapted from ASD testing code

%% Setup
% Get data 
X = imread("datasets/lena.png");
M=double(X)/255;
r = 50;
[U, S, V] = svd(M);
% clear X;
% M = U*S*V';

% m is number of rows, n is number of columns, delta is percent of entries
% the algorithm sees, p is number of observations
[m, n] = size(M);
delta = 0.3;
p = round(delta*m*n);

% At first, A contains p 1s followed by m * n - p zeros
% Then A gets randomly shuffled, and finally turned into an m * n array
A = zeros(m*n,1);
A(1:p) = ones(p,1);
[~,ind] = sort(randn(m*n,1));
A = A(ind);
A = reshape(A,m,n);

% Generate rows and column indices (I and J), turn them into linear indices
% (Omega), then create an obvservations array which is just the values the
% algo can see
[I, J] = find(A);
Omega = sub2ind([m n], I, J);
obvservations=M(Omega);

%% Run algorithm
max_time = 60;
output_file = fopen('output/lena.txt', 'w');
alg_names = ["R3MC", "ScaledASD", "ASD", "NIHT_Matrix", "CGIHT_Matrix", "ScaledGD", "LMaFit"];
run_test(A, obvservations, r, alg_names, max_time, output_file, M, ones(size(M)));
fclose(output_file);
