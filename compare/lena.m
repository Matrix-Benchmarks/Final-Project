%% Intro
% Adapted from ASD testing code
% TODO: Get other algorithms besides ASD working

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

% Display what the algorithm will see, which is just M masked to A (<M, A>)
% Note this is basically equivalent to data above, just in a 2d array
% format for easy visualization
% figure(1)
% imshow(A.*M)
% drawnow

%% Run algorithm

alg_names = {'ScaledASD', 'ASD'};
[Xr,outs,alg_names] = run_MC_algos(A,obvservations,r,alg_names,default_params(cond_nr, m, n));

%% Print frobenius info to file
output_file = fopen('output/lena.txt','w');
for algo_num = 1:size(Xr, 2)
    for iterate = 1:size(Xr{algo_num}, 2)
       current = Xr{algo_num}{iterate};
       time = outs{algo_num}.time(iterate);
       try 
           matrix_iterate =  current{1} * transpose(current{2});
       catch
           matrix_iterate = current.U * transpose(current.V);
       end
       fprintf(output_file, "%s %i %d %d\n", alg_names{algo_num}, iterate, time(1), norm(matrix_iterate - M, 'fro'));
    end
end