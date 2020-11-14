output_file = fopen('output/test.txt', 'w');
% TODO: Fix algos that aren't quite finishing
% for matrix_size = [200, 2000, 20000, 200000]
max_time = 30;
for matrix_size = [10000]
    for rank = [2, 20, 200]
%         for fraction_shown = [0.05, 0.1, 0.2, 0.4]
        for fraction_shown = [0.1, 0.2, 0.4]
            for condition_number = [1.5]
                if (matrix_size^2 * fraction_shown > (2 * matrix_size - rank) * rank) % Information theoretic limit per https://arxiv.org/abs/1504.04970
                    fprintf(output_file, "Parameters: %i %i %i %f\n", matrix_size, rank, fraction_shown * 100, condition_number);
                    run_random(rank, matrix_size, fraction_shown, condition_number, output_file, max_time)
                    clearvars -except output_file matrix_size rank fraction_shown max_time
                end
            end
        end
        
    end
end
fclose(output_file);