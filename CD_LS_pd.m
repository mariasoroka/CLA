function [x, log_resid, log_x, log_time, log_flops, n_iter] = CD_LS_pd(A, b, e, max_time, verbose,mysol)
    [m, n] = size(A);

    [log_resid, log_x, log_time, log_flops] = deal([]);

    x = zeros(n, 1);
    n_iter = 0;
    time = 0;
    residue = (dot(x-mysol,A*(x-mysol))/dot(mysol,A*mysol))^0.5;
    I = eye(n); 
 
    prob_array = diag(A) / trace(A);
    if issparse(A) == true
        prob_array = full(prob_array);
    end
    dpdf = makedist('Multinomial', 'Probabilities', prob_array);

    flops(0);
    tic;
    
    while residue > e && time < max_time
        n_iter = n_iter + 1;
        n_col = random(dpdf);
        %disp(A(n_col,:))
        x = x - (dot(A(n_col, :),x)-b(n_col)) / A(n_col, n_col).*I(:, n_col);

        residue = (dot(x-mysol,A*(x-mysol))/dot(mysol,A*mysol))^0.5;

        if verbose == true
            if issparse(A) == true
                addflops(flops_spmul(A, x));
            else
                addflops(flops_mul(A, x));
            end

            addflops(n); % for the dot product
            addflops(n); % for the substraction
            addflops(n*flops_div); % for the division in each term
            addflops(2*n+2);  % to do this operation (norm(A(i,:))^2)
            addflops(n); % for the scalar multiplication for myvec

            log_resid = [log_resid, residue];
            log_x = [log_x, x];
            log_time = [log_time, toc];
            log_flops = [log_flops, flops];
        end
        time = toc;
    end
end