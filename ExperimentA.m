%% Generate consistent system
m = 1000;
n = 500;
A = rand(m, n);
x_gt = rand(n, 1);
b = A * x_gt;
%% Solve system using Kaczmarz
[x, log_resid, log_x, n_iter] = Kaczmarz(A, b, 10e-3, m * 100, true);
%% Get ||x - x_gt||
resid = vecnorm(log_x - x_gt);
%% Plot residuals
iters = 1:n_iter;
plot(iters, log10(resid(1:n_iter)))
hold on
%% Solve system using simple randomized Kaczmarz
[x, log_resid, log_x, n_iter] = KaczmarzRandom(A, b, 10e-3, m * 100, true, true);
%% Get ||x - x_gt||
resid = vecnorm(log_x - x_gt);
%% Plot residual
iters = 1:n_iter;
plot(iters, log10(resid(1:n_iter)))
hold on
%% Solve system using simple randomized Kaczmarz
[x, log_resid, log_x, n_iter] = KaczmarzRandom(A, b, 10e-3, m * 100, true, false);
%% Get ||x - x_gt||
resid = vecnorm(log_x - x_gt);
%% Plot residual
iters = 1:n_iter;
plot(iters, log10(resid(1:n_iter)))
%% Labeling
xlabel("number of iterations")
ylabel("log(||x_k - x||)")
legend("Kaczmarz method", "Simple randomized Kaczmarz", "Randomized Kaczmarz")
hold off
%%