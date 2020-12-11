%% Get Dataset
T = readtable('datasets/ratings.dat');
M = zeros(4000);
[num_rows, ~] = size(T);
total_observations = 0;
last_user = -1;
for i = 1:num_rows
   user = T{i, 1};
   movie = T{i, 3};
   rating = T{i, 5}; 
   if mod(user, 40) == 0 && user ~= last_user
       disp(user / 40);
       last_user = user;
   end
   if user < 4000 && movie < 4000
      M(user, movie) = rating;
      total_observations = total_observations + 1;
   end
end

%% Prepare for running algos
r = 10;
positions = java.util.HashMap;
percent_shown = 0.3;
p = total_observations * percent_shown;

while positions.size() < p
    test = randi([1, total_observations]);
    while positions.containsKey(test) && M(test) ~= 99
        test = randi([1, total_observations]);
    end
    positions.put(test, 0);
end

filter = M;
filter(filter~=0) = 1;

array = positions.keySet().toArray();
A = zeros(size(M));
for i = 1:p
    obs_index = array(i);
    user = T{obs_index, 1};
    movie = T{obs_index, 3};
    rating = T{obs_index, 5}; 
    A(user, movie) = rating;
end

% Generate rows and column indices (I and J), turn them into linear indices
% (Omega), then create an obvservations array which is just the values the
% algo can see
[I, J] = find(A);
Omega = sub2ind([4000 4000], I, J);
observations=M(Omega);

%% Run algos
max_time = 60;
output_file = fopen('output/results/movielens.txt', 'a');
alg_names = ["MatrixIRLS"];
run_test(A, observations, r, alg_names, max_time, output_file, M, filter, sparse(A));
fclose(output_file);
