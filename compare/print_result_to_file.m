function print_result_to_file(algo_name, time_info, iteration_info, output_file, goal_matrix) 
    for iterate = 1:size(iteration_info, 2)
        current = iteration_info{iterate};
        time = time_info.time(iterate);

        try
            matrix_iterate =  current{1} * transpose(current{2});
        catch
            % Get dense matrix from sparse representation
            matrix_iterate = get_densemat_from_compact(current, A); 
            disp(norm(matrix_iterate - goal_matrix, 'fro'));
        end

        try
            split_name = split(algo_name);
            fprintf(output_file, "%s %i %d %d\n", split_name{1}, iterate, time(1), norm(matrix_iterate - goal_matrix, 'fro'));
        catch
                  
                fprintf("Printing to file failed for algorithm number %i iteration %i\n", algo_num, iterate)
        end
    end
