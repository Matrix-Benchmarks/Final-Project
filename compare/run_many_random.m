output_file = fopen('output/results/randoms.txt', 'w');
rng(42)
% TODO: Fix algos that aren't quite finishing
max_time = 60;
for matrix_size = [1000, 10000]
    for rank = [2, 20, 200]
        for fraction_shown = [0.05, 0.1, 0.2, 0.4]
            for condition_number = [2, 200, 20000]
                if (matrix_size^2 * fraction_shown > (2 * matrix_size - rank) * rank) % Information theoretic limit per https://arxiv.org/abs/1504.04970
                    fprintf(output_file, "Parameters: %i %i %i %i\n", matrix_size, rank, fraction_shown * 100, condition_number);
                    run_random(rank, matrix_size, fraction_shown, condition_number, output_file, max_time)
                    clearvars -except output_file matrix_size rank fraction_shown max_time
                end
            end
        end
        
    end
end
fclose(output_file);