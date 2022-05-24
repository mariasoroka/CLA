e = 10^(-2);
t = 100;
k = 10;

%%
a = 2;
mymean = a;
for j=1:k
     mymean = (1/(j+1))*(j*mymean+a);
end
disp(mymean)
%% Dense case
n = 100;
A = hilb(n);
cond(A)

% Computing ground truth
x_gt = randn(n, 1);
b = A * x_gt;


% Solve system using Gauss_LS_pd
[x_G_LS_dense, log_resid_G_LS_dense, log_x_G_LS_dense, log_time_G_LS_dense, log_flops_G_LS_dense, n_Gausspd] = Gauss_pd(A, b, e, t, true);
% Solve system using CD_LS_pd
[x_CD_LS_dense, log_resid_CD_LS_dense, log_x_CD_LS_dense, log_time_CD_LS_dense, log_flops_CD_LS_dense, n_LSpd] = CD_LS_pd(A, b, e, t, true);
mynGauss = n_Gausspd;
mynCD = n_LSpd;

mean_log_resid_G_LS_dense = log_resid_G_LS_dense;
mean_log_time_G_LS_dense  = log_time_G_LS_dense;
mean_log_flops_G_LS_dense = log_flops_G_LS_dense;

mean_log_resid_CD_LS_dense = log_resid_CD_LS_dense;
mean_log_time_CD_LS_dense = log_time_CD_LS_dense;
mean_log_flops_CD_LS_dense= log_flops_CD_LS_dense;

disp(mean_log_time_CD_LS_dense)
disp(mean_log_time_G_LS_dense)
%%
for j=1:k
    % Solve system using Gauss_LS_pd
    [x_G_LS_dense, log_resid_G_LS_dense, log_x_G_LS_dense, log_time_G_LS_dense, log_flops_G_LS_dense, n_Gausspd] = Gauss_pd(A, b, e, t, true);
    % Solve system using CD_LS_pd
    [x_CD_LS_dense, log_resid_CD_LS_dense, log_x_CD_LS_dense, log_time_CD_LS_dense, log_flops_CD_LS_dense, n_LSpd] = CD_LS_pd(A, b, e, t, true);
    mynGauss = min(n_Gausspd,mynGauss);
    disp(mynGauss) 
    disp(log_time_G_LS_dense(end))
    %disp(log_time_G_LS_dense(mynGauss))
    %disp(mean_log_time_G_LS_dense(mynGauss))
    %disp(log_resid_G_LS_dense)
    %disp(log_time_G_LS_dense)
    mean_log_resid_G_LS_dense =    (1/(j+1))*(j*mean_log_resid_G_LS_dense(1:mynGauss)+log_resid_G_LS_dense(1:mynGauss));
    mean_log_time_G_LS_dense  = (1/(j+1))*(j*mean_log_time_G_LS_dense(1:mynGauss)+log_time_G_LS_dense(1:mynGauss));
    mean_log_flops_G_LS_dense = (1/(j+1))*(j*mean_log_flops_G_LS_dense(1:mynGauss)+log_flops_G_LS_dense(1:mynGauss));
    
    %disp("changealgo")
    %disp(mynGauss) 
    %disp(log_time_CD_LS_dense(mynGauss))
    %disp(mean_log_time_CD_LS_dense(mynGauss))
    mynCD = min(mynCD,n_LSpd);
    mean_log_resid_CD_LS_dense = (1/(j+1))*(j*mean_log_resid_CD_LS_dense(1:mynCD)+log_resid_CD_LS_dense(1:mynCD));
    mean_log_time_CD_LS_dense = (1/(j+1))*(j*mean_log_time_CD_LS_dense(1:mynCD)+log_time_CD_LS_dense(1:mynCD));
    mean_log_flops_CD_LS_dense= (1/(j+1))*(j*mean_log_flops_CD_LS_dense(1:mynCD)+log_flops_CD_LS_dense(1:mynCD));
    
end






%% Sparse case
n = 1000;
density = 1/log(n^2);
rc = 1/n;
A = sprandsym(n, density, rc, 1);
%% Computing ground truth
x_gt = randn(n, 1);
b = A * x_gt;
%% Solve system using Gauss_LS
[x_G_LS_sparse, log_resid_G_LS_sparse, log_x_G_LS_sparse, log_time_G_LS_sparse, log_flops_G_LS_sparse, ~] = Gauss_pd(A, b, e, t, true);
%% Solve system using CD_LS
[x_CD_LS_sparse, log_resid_CD_LS_sparse, log_x_CD_LS_sparse, log_time_CD_LS_sparse, log_flops_CD_LS_sparse, ~] = CD_LS(A, b, e, t, true);

%%
% Plot time - dense
subplot(2, 2, 1)
semilogy(mean_log_time_G_LS_dense, mean_log_resid_G_LS_dense);
hold on
semilogy(mean_log_time_CD_LS_dense, mean_log_resid_CD_LS_dense);
xlabel("time (s)")
ylabel("error")
legend("Gauss LS pd", "CD LS pd")
hold off
%%
% Plot flops - dense
subplot(2, 2, 2)
semilogy(mean_log_flops_G_LS_dense, mean_log_resid_G_LS_dense);
hold on
semilogy(mean_log_flops_CD_LS_dense, mean_log_resid_CD_LS_dense);
xlabel("flops")
ylabel("error")
legend("Gauss LS pd", "CD LS pd")
hold off

%% Plot time - sparse
subplot(2, 2, 3)
semilogy(log_time_G_LS_sparse, log_resid_G_LS_sparse);
hold on
semilogy(log_time_CD_LS_sparse, log_resid_CD_LS_sparse);
xlabel("time")
ylabel("error")
legend("Gauss LS pd", "CD LS pd")
hold off
%% Plot flops - sparse
subplot(2, 2, 4)
semilogy(log_flops_G_LS_sparse, log_resid_G_LS_sparse);
hold on
semilogy(log_flops_CD_LS_sparse, log_resid_CD_LS_sparse);
xlabel("flops")
ylabel("error")
legend("Gauss LS pd", "CD LS pd")
hold off