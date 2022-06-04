%% Generate consistent system
m = 1000;
n = 500;
A = randn(m, n);
% for i=1:m
%     if(mod(i, 10) == 0)
%         A(i, :) = A(i, :) * 100;
%     end
% end
x_gt = randn(n, 1);
b = A * x_gt;
%% Solve system using Kaczmarz
[x, log_resid, log_x, n_iter] = Kaczmarz(A, b, 10e-3, 20 * m, true);
%% Plot residuals
iters = 1:n_iter;
semilogy(iters, log_resid(1:n_iter))
hold on
%% Solve system using simple randomized Kaczmarz
[x, log_resid, log_x, n_iter] = KaczmarzRandom(A, b, 10e-3, 20 * m, true, true);
%% Plot residual
iters = 1:n_iter;
semilogy(iters, log_resid(1:n_iter))
hold on
%% Solve system using proper randomized Kaczmarz
[x, log_resid, log_x, n_iter] = KaczmarzRandom(A, b, 10e-3, 20 * m, true, false);
%% Plot residual
iters = 1:n_iter;
semilogy(iters, log_resid(1:n_iter))
%% Labeling
xlabel("number of iterations")
ylabel("log(||Ax - b||)")
legend("Kaczmarz method", "Simple randomized Kaczmarz", "Randomized Kaczmarz")
hold off