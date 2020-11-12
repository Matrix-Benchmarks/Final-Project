% TODO: Add condition number
% TODO: Not square matrices?
% TODO: Give wrong rank to algos on purpose?
% TODO: Noise?
% This function prints algo_name, iteration number, time, frob distance to
% the file pointer output_file for each algo and iteration
function run_random(rank, matrix_size, fraction_shown, output_file, max_time)
    

    %% Setup random matrix 
    goal_matrix = randn(matrix_size);
    [U, S, V] = svds(goal_matrix, rank);
    goal_matrix = U * S * V';

    
    % Generate condition number, needed for some algos
    singular_vals = find(goal_matrix);
    cond_nr = singular_vals(rank) / singular_vals(1);
    
    num_observations = round(fraction_shown*matrix_size*matrix_size);
    
    % At first, the mask contains num_observations 1s followed by size * size - num_observations zeros
    % Then the mask gets randomly shuffled, and finally turned into a size * size array
    mask = zeros(matrix_size*matrix_size,1);
    mask(1:num_observations) = ones(num_observations,1);
    [~,ind] = sort(randn(matrix_size*matrix_size,1));
    mask = mask(ind);
    mask = reshape(mask,matrix_size,matrix_size);
    
    % Generate rows and column indices (I and J), turn them into linear indices
    % (Omega), then create an obvservations array which is just the values the
    % algo can see
    [I, J] = find(mask);
    Omega = sub2ind([matrix_size matrix_size], I, J);
    obvservations = goal_matrix(Omega);
    

    %% Run matrix through algorithms
    alg_names = ["R3MC", "ScaledASD", "ASD", "NIHT_Matrix", "CGIHT_Matrix", "ScaledGD", "LMaFit"];
    [iteration_info, time_info] = run_test(mask, obvservations, rank, alg_names, max_time);
        
    %% Print frobenius info to file
    for algo_num = 1:size(iteration_info, 2)
        for iterate = 1:size(iteration_info{algo_num}, 2)
            current = iteration_info{algo_num}{iterate};
            time = time_info{algo_num}.time(iterate);
            try
                matrix_iterate =  current{1} * transpose(current{2});
            catch
%                 disp(current)
                % Get dense matrix from sparse representation
                matrix_iterate = get_densemat_from_compact(current, A); 
                disp(norm(matrix_iterate - goal_matrix, 'fro'));
            end

            try
                split_name = split(alg_names{algo_num});
                fprintf(output_file, "%s %i %d %d\n", split_name{1}, iterate, time(1), norm(matrix_iterate - goal_matrix, 'fro'));
            catch
%                 fprintf("Printing to file failed for algorithm number %i iteration %i\n", algo_num, iterate)
            end
        end
    end
end