output_file = fopen('output/test.txt', 'a');
max_time = 30;
for matrix_size = [200, 2000, 20000, 200000]
    for rank = [2, 20, 200]
        for fraction_shown = [0.05, 0.1, 0.2, 0.4, 0.8]
            run_random(rank, matrix_size, fraction_shown, output_file, max_time)
        end
    end
end
fclose(output_file);