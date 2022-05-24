function [x, log_resid, log_x, n_iter] = Kaczmarz(A, b, e, max_n_iter, verbose)
    [m, n] = size(A);
    x = randn(n, 1);
    n_iter = 0;
    row_n = 1;
    residue = norm(A * x - b);
    log_resid = zeros(1, max_n_iter);
    log_x = zeros(n, max_n_iter);
    while residue > e && n_iter < max_n_iter
        n_iter = n_iter + 1;
        if verbose == true
            log_resid(n_iter) = residue;
            log_x(:, n_iter) = x;
        end
        x = x + (b(row_n) - A(row_n, :) * x) / (norm(A(row_n, :))^2) * A(row_n, :)';
        row_n = row_n + 1;
        if row_n > m 
            row_n = 1;
        end
        residue = norm(A * x - b);
    end
end