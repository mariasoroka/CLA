% Gaussian Kaczmarz (A not necessarily positive definite)
%theoretical part: B=I; sigma = I
function [x, log_resid, log_x, log_time, log_flops, n_iter] = GK(A, b, e, max_time, verbose)
    [m, n] = size(A);

    [log_resid, log_x, log_time, log_flops] = deal([]);

    x = rand(n, 1);
    n_iter = 0;
    time = 0;
    residue = norm(A * x - b) / norm(b);

    id = eye(m);
    mu = zeros(m,1);

    flops(0);
    tic;
    
    while residue > e && time < max_time
        n_iter = n_iter + 1;
        eta = mvnrnd(mu, id);
        eta = eta';
        x = x - eta' * (A*x - b) / (norm(A'*eta)^2) .* A' * eta;

        residue = norm(A * x - b) / norm(b);

        if verbose == true
            if issparse(A) == true
                addflops(2*flops_spmul(A', eta));
                addflops(flops_spmul(A, x));
            else
                addflops(2*flops_mul(A', eta));
                addflops(flops_mul(A, x));
            end

            addflops(n^2); %to compute the transpose
            addflops(n); % for the substraction
            addflops(n); % for the vector product
            addflops(n * flops_div); % for the division (on each component)
            addflops(2*n + 2);  % to compute the norm squared
            addflops(n); % for the scalar multiplication for myvec
            log_resid = [log_resid, residue];
            log_x = [log_x, x];
            log_time = [log_time, toc];
            log_flops = [log_flops, flops];
        end
        time = toc;
    end
end