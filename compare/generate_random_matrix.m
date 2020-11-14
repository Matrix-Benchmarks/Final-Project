function matrix=generate_random_matrix(type, matrix_size, rank, cond_num)

% Some code taken from sample_X0_lowrank.m, written by Christian Kuemmerle

if (type == "LOG")
    fct = @(l) cond_num.*exp(-log(cond_num).*([1:l]-1)./(rank-1));
    S = fct(rank);
%     S = cond_num.*exp(-log(cond_num).*([1:rank]-1)./(rank-1));
else
    disp("Unsupported singular value distribution!")
    quit(1)
end

U = orth(randn(matrix_size, rank)); % Random orthonormal n * r
V = orth(randn(matrix_size, rank))'; % Random orthonormal r * n
matrix = (U.*S) * V;
end