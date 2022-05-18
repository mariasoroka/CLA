e = 10^(-4);
%% Dense case
n = 100;
A = hilb(n);
cond(A)
%% Computing ground truth
x_gt = rand(n, 1);
b = A * x_gt;
%% Solve system using Gauss_LS_pd
[x_G_LS_dense, log_resid_G_LS_dense, log_x_G_LS_dense, log_time_G_LS_dense, log_flops_G_LS_dense, ~] = Gauss_pd(A, b, e, 300, true);
%% Solve system using CD_LS_pd
[x_CD_LS_dense, log_resid_CD_LS_dense, log_x_CD_LS_dense, log_time_CD_LS_dense, log_flops_CD_LS_dense, ~] = CD_LS_pd(A, b, e, 300, true);




%% Sparse case
n = 1000;
density = 1/log(n^2);
rc = 1/n;
A = sprandsym(n, density, rc, 1);
%% Computing ground truth
x_gt = rand(n, 1);
b = A * x_gt;
%% Solve system using Gauss_LS
[x_G_LS_sparse, log_resid_G_LS_sparse, log_x_G_LS_sparse, log_time_G_LS_sparse, log_flops_G_LS_sparse, ~] = Gauss_LS(A, b, e, 300, true);
%% Solve system using CD_LS
[x_CD_LS_sparse, log_resid_CD_LS_sparse, log_x_CD_LS_sparse, log_time_CD_LS_sparse, log_flops_CD_LS_sparse, ~] = CD_LS(A, b, e, 300, true);



%% Plot time - dense
subplot(2, 2, 1)
plot(log_time_G_LS_dense, log_resid_G_LS_dense);
hold on
plot(log_time_CD_LS_dense, log_resid_CD_LS_dense);
xlabel("time")
ylabel("error")
legend("Gauss LS pd", "CD LS pd")
hold off
%% Plot flops - dense
subplot(2, 2, 2)
plot(log_flops_G_LS_dense, log_resid_G_LS_dense);
hold on
plot(log_flops_CD_LS_dense, log_resid_CD_LS_dense);
xlabel("flops")
ylabel("error")
legend("Gauss LS pd", "CD LS pd")
hold off
%% Plot time - sparse
subplot(2, 2, 3)
plot(log_time_G_LS_sparse, log_resid_G_LS_sparse);
hold on
plot(log_time_CD_LS_sparse, log_resid_CD_LS_sparse);
xlabel("time")
ylabel("error")
legend("Gauss LS pd", "CD LS pd")
hold off
%% Plot flops - sparse
subplot(2, 2, 4)
plot(log_flops_G_LS_sparse, log_resid_G_LS_sparse);
hold on
plot(log_flops_CD_LS_sparse, log_resid_CD_LS_sparse);
xlabel("flops")
ylabel("error")
legend("Gauss LS pd", "CD LS pd")
hold off