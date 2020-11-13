output_file = fopen('output/test.txt', 'a');
% TODO: Fix algos that aren't quite finishing
% for matrix_size = [200, 2000, 20000, 200000]
for matrix_size = [20000]
    for rank = [2, 20, 200]
%         for fraction_shown = [0.05, 0.1, 0.2, 0.4]
        for fraction_shown = [0.1, 0.2, 0.4]
            if (matrix_size^2 * fraction_shown > (2 * matrix_size - rank) * rank) % Information theoretic limit per https://arxiv.org/abs/1504.04970
                max_time = 30;
                fprintf(output_file, "Parameters: %i %i %i\n", matrix_size, rank, fraction_shown * 100);
                run_random(rank, matrix_size, fraction_shown, output_file, max_time)
                clearvars -except output_file matrix_size rank fraction_shown
            end
        end
    end
end
fclose(output_file);