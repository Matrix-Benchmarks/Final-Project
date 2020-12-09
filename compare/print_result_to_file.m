function print_result_to_file(algo_name, time_info, iteration_info, output_file, goal_matrix, filter, omega) 
    for iterate = 1:size(iteration_info, 2)
        current = iteration_info{iterate};
        time = time_info.time(iterate);

        try
            matrix_iterate =  current{1} * transpose(current{2});
        catch
            try 
            matrix_iterate = get_densemat_from_compact(current, omega); 
            catch
                fprintf("Parsing failed for algorithm name %s iteration %i\n", algo_name, iterate);
                continue;
            end
        end
        
        try
            split_name = split(algo_name);
            fprintf(output_file, "%s %i %d %d\n", split_name{1}, iterate, time(1), norm(dot(matrix_iterate, filter) - dot(goal_matrix, filter), 'fro'));
        catch
            fprintf("Printing to file failed for algorithm name %s iteration %i\n", algo_name, iterate)
        end
    end
