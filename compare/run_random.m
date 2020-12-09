% TODO: Add condition number
% TODO: Not square matrices?
% TODO: Give wrong rank to algos on purpose?
% TODO: Noise?
% This function prints algo_name, iteration number, time, frob distance to
% the file pointer output_file for each algo and iteration
function run_random(rank, matrix_size, fraction_shown, condition_number, output_file, max_time)
    

    %% Setup random matrix 
    goal_matrix = generate_random_matrix("LOG", matrix_size, rank, condition_number);
    fprintf("Matrix built, rank %i, fraction shown %f, condition_number %i, %i max_time\n", rank, fraction_shown, condition_number, max_time)

    num_observations = round(fraction_shown*matrix_size*matrix_size);
    % At first, the mask contains num_observations 1s followed by size * size - num_observations zeros
    % Then the mask gets randomly shuffled, and finally turned into a size * size array
    mask = zeros(matrix_size*matrix_size,1);
    fprintf("Mask Made\n")
    one_positions =  randperm(matrix_size * matrix_size, num_observations);
    fprintf("Positions generated\n")
    mask(one_positions) = 1;
    fprintf("Mask Done\n")
    mask = reshape(mask,matrix_size,matrix_size);
    
    % Generate rows and column indices (I and J), turn them into linear indices
    % (Omega), then create an obvservations array which is just the values the
    % algo can see
    [I, J] = find(mask);
    Omega = sub2ind([matrix_size matrix_size], I, J);
    obvservations = goal_matrix(Omega);
    
    %% Run matrix through algorithms
    disp("Running algos")
    alg_names = ["MatrixIRLS"];
    run_test(mask, obvservations, rank, alg_names, max_time, output_file, goal_matrix, ones(size(goal_matrix)), sparse(mask));
        
end