m = 500;
n = 50; 
e = 10^(-4);

%% Dense case
A = randn(m, n);
%% Computing ground truth
x_gt = rand(n, 1);
b = A * x_gt;
k = cond(A)
%% Compute theoretical convergence rate
rho = 1 - eigs(A.' * A, 1, "smallestabs") / norm(A, "fro")^2;
%% Solve system using Gauss_LS
[x_G_LS_dense, log_resid_G_LS_dense, log_x_G_LS_dense, log_time_G_LS_dense, log_flops_G_LS_dense, n_iters_G_LS_dense] = Gauss_LS(A, b, e, 300, true);
%% Compute quantiles
n_runs = 100;
max_n_iters = 2000;
res = zeros(n_runs, max_n_iters);
num_iters = zeros(1, n_runs);
for i = 1:n_runs
    [x, log_resid, log_x, log_time, log_flops, n_iters] = Gauss_LS(A, b, e, 15, true);
    num_iters(i) = n_iters;
    res(i, 1:min(max_n_iters, n_iters)) = log_resid(1:min(max_n_iters, n_iters));
end
all_iter = min(num_iters);
Q = quantile(res(:, 1:all_iter), [0.05, 0.95], 1);
%% Plot residuals - time
iters = 1:n_iters_G_LS_dense;
all_iters = 1:all_iter;
plot(iters, log10(log_resid_G_LS_dense));
hold on
plot(iters, iters * log10(rho) + log10(log_resid_G_LS_dense(1)));
hold on
plot(all_iters, log10(Q(1, :)));
hold on
plot(all_iters, log10(Q(2, :)));
xlabel("iterations")
ylabel("error")
legend("Gauss LS", "theoretical")
hold off