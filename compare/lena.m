%% Intro
% Adapted from ASD testing code

%% Setup
% Get data and project to closest 50 rank matrix M
X = imread("datasets/lena.png");
X=double(X)/255;
r = 50;
[U, S, V] = svds(X,r);
clear X;
M = U*S*V';

% Generate condition number, needed for some algos
singular_vals = find(S);
cond_nr = singular_vals(50) / singular_vals(1);


% m is number of rows, n is number of columns, delta is percent of entries
% the algorithm sees, p is number of observations
[m, n] = size(M);
delta = 0.35;
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

% Work: "R3MC", "ScaledASD", "ASD", "NIHT_Matrix", "CGIHT_Matrix"
% TODO: LRGeomCG, MatrixIRLS, other IRLS?, ScaledGD, other Reimann ones?
alg_names = ["LMaFit", "R3MC", "ScaledASD", "ASD", "NIHT_Matrix", "CGIHT_Matrix"];
[Xr,outs] = run_test(A,obvservations,r,alg_names,60);

%% Print frobenius info to file
output_file = fopen('output/lena.txt','w');
for algo_num = 1:size(Xr, 2)
    for iterate = 1:size(Xr{algo_num}, 2)
        current = Xr{algo_num}{iterate};
        time = outs{algo_num}.time(iterate);
        try
            matrix_iterate =  current{1} * transpose(current{2});
        catch
            % Get dense matrix from sparse representation
            matrix_iterate = get_densemat_from_compact(current, A); 
            disp(norm(matrix_iterate - M, 'fro'));
        end
        
        try
            split_name = split(alg_names{algo_num});
            fprintf(output_file, "%s %i %d %d\n", split_name{1}, iterate, time(1), norm(matrix_iterate - M, 'fro'));
        catch
        end
    end
end