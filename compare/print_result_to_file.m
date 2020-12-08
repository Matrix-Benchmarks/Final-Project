function print_result_to_file(algo_name, time_info, iteration_info, output_file, goal_matrix, filter) 
    for iterate = 1:size(iteration_info, 2)
        current = iteration_info{iterate};
        time = time_info.time(iterate);

%         try
%             x = randi([1 size(goal_matrix, 1)], 1, sample_size);
%             y = randi([1 size(goal_matrix, 2)], 1, sample_size);
%             sparse_test_values = partXY(current{1}', current{2}', x, y, sample_size).';
%             linear_inds = sub2ind(size(goal_matrix), x, y);
%             sparse_goal_values = goal_matrix(linear_inds);
%         catch
%             % Get dense matrix from sparse representation
% %             matrix_iterate = get_densemat_from_compact(current, A); 
% %             disp(norm(matrix_iterate - goal_matrix, 'fro'));
%         end

        try
            matrix_iterate =  current{1} * transpose(current{2});
            split_name = split(algo_name);
            fprintf(output_file, "%s %i %d %d\n", split_name{1}, iterate, time(1), norm(dot(matrix_iterate, filter) - dot(goal_matrix, filter), 'fro'));
        catch
                  
                fprintf("Printing to file failed for algorithm name %s iteration %i\n", algo_name, iterate)
        end
    end
