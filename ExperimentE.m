e = 10^(-4);
t = 100;
%% Dense case
n = 100;
A = hilb(n);
cond(A)
%% Computing ground truth
x_gt = rand(n, 1);
b = A * x_gt;
%% Solve system using Gauss_LS_pd
[x_G_LS_dense, log_resid_G_LS_dense, log_x_G_LS_dense, log_time_G_LS_dense, log_flops_G_LS_dense, ~] = Gauss_pd(A, b, e, t, true,x_gt);
%% Solve system using CD_LS_pd
[x_CD_LS_dense, log_resid_CD_LS_dense, log_x_CD_LS_dense, log_time_CD_LS_dense, log_flops_CD_LS_dense, ~] = CD_LS_pd(A, b, e, t, true,x_gt);




%% Sparse case
n = 1000;
density = 1/log(n^2);
rc = 1/n;
A = sprandsym(n, density, rc, 1);
%% Computing ground truth
x_gt = randn(n, 1);
b = A * x_gt;
%% Solve system using Gauss_LS
[x_G_LS_sparse, log_resid_G_LS_sparse, log_x_G_LS_sparse, log_time_G_LS_sparse, log_flops_G_LS_sparse, ~] = Gauss_pd(A, b, e, t, true,x_gt);
%% Solve system using CD_LS
[x_CD_LS_sparse, log_resid_CD_LS_sparse, log_x_CD_LS_sparse, log_time_CD_LS_sparse, log_flops_CD_LS_sparse, ~] = CD_LS_pd(A, b, e, t, true,x_gt);



%% Plot time - dense
subplot(2, 2, 1)
semilogy(log_time_G_LS_dense, (100/log_resid_G_LS_dense(1))*log_resid_G_LS_dense);
hold on
semilogy(log_time_CD_LS_dense,(100/log_resid_CD_LS_dense(1))*log_resid_CD_LS_dense);
xlabel("time")
ylabel("error")
legend("Gauss LS pd", "CD LS pd")
hold off
%% Plot flops - dense
subplot(2, 2, 2)
semilogy(log_flops_G_LS_dense, (100/log_resid_G_LS_dense(1))*log_resid_G_LS_dense);
hold on
semilogy(log_flops_CD_LS_dense, (100/log_resid_CD_LS_dense(1))*log_resid_CD_LS_dense);
xlabel("flops")
ylabel("error")
legend("Gauss LS pd", "CD LS pd")
hold off
%% Plot time - sparse
subplot(2, 2, 3)
semilogy(log_time_G_LS_sparse, (100/log_resid_G_LS_sparse(1))*log_resid_G_LS_sparse);
hold on
semilogy(log_time_CD_LS_sparse, (100/log_resid_CD_LS_sparse(1))*log_resid_CD_LS_sparse);
xlabel("time")
ylabel("error")
legend("Gauss LS pd", "CD LS pd")
hold off
%% Plot flops - sparse
subplot(2, 2, 4)
semilogy(log_flops_G_LS_sparse, (100/log_resid_G_LS_sparse(1))*log_resid_G_LS_sparse);
hold on
semilogy(log_flops_CD_LS_sparse, (100/log_resid_CD_LS_sparse(1))*log_resid_CD_LS_sparse);
xlabel("flops")
ylabel("error")
legend("Gauss LS pd", "CD LS pd")
hold off