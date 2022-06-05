function [x, log_resid, log_x, n_iter] = KaczmarzRandom(A, b, e, max_n_iter, verbose, simple)
    [m, n] = size(A);
    frob_norm = norm(A,"fro"); 
    if simple == false
        prob_array = vecnorm(A.').^2 ./ frob_norm^2;
    else
        prob_array = ones(1, m)./m;
    end
    prob_array = full(prob_array);

    dpdf = makedist('Multinomial', 'Probabilities', prob_array);
    x = zeros(n, 1);
    n_iter = 0;
    residue = norm(A*x-b);
    log_resid = zeros(1, max_n_iter);
    log_x = zeros(n, max_n_iter);
    while residue > e && n_iter < max_n_iter
        n_iter = n_iter + 1;
        row_n = random(dpdf);
        if verbose == true
            log_resid(n_iter) = residue;
            log_x(:, n_iter) = x;
        end
        x = x + (b(row_n) - A(row_n, :) * x) / (norm(A(row_n, :))^2) * A(row_n, :)';
        residue = norm(A * x - b);
    end
end

