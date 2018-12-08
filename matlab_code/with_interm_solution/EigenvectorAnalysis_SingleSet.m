%% number of cells in the model
n_cells = 3;

%% initialize m_n, m_c, p_c and p_n for each cell as well as p_ECF and V_ECF
index = 1:4:n_cells * 4;
all_values = zeros(1, n_cells * 4 + 1);
all_values(index) = 10.9;
all_values(index + 1) = 1.88;
all_values(index + 2) = 129.21;
all_values(index + 3) = 3531.75;

% initialize p_ECF (concentration of PER proteins in Extracellular Fuild) 
all_values(end) = 18;

% initialize the volumn of Extracellular Fuild
V_ECF = 1e-4;

%% initialize K, K_ex, K_im, alpha, beta, r_ex, r_im
%  same for each cell
K = 200;
alpha = 0.04;
beta = 0.04;
lambda_ECF = 0.46;
V_n = 1e-5;
V_c = 5e-5;
%% initialize gamma_m = delta_m = delta_p = gamma_p = lambda as well as r
%  different for each cell
r = zeros(1, n_cells);
r(1) = 10;
r(2) = 10;
r(3:end) = 10;
lambda = zeros(1, n_cells);
lambda(1) = 1.2 * pi / 220;
lambda(2) = 0.8 * pi / 110;
lambda(3:end) = 1 * pi / 110;
%%

for j = 1:2
    if j == 2
        temp = lambda_ECF;
        lambda_ECF = 0;
    end
    fun = @(x)SolveSteady(x, n_cells, alpha, V_n, K, r, lambda, V_c, beta, lambda_ECF, V_ECF);
    
    x0 = zeros(1, n_cells * 4 + 1);
    options = optimoptions('fsolve'); 
    options.MaxIterations = 100000;
    options.MaxFunctionEvaluations = 100000000;
    x = fsolve(fun, x0, options);

    a = zeros(1, n_cells);
    a_ex = zeros(1, n_cells);

    disp(x);

    for i = 1:n_cells
        a(i) = - (alpha / V_n) * r(i) * (K / (K + x((i - 1)*4 + 4))) .^ r(i) * (1 / (K + x((i - 1)*4 + 4)));
    end

    V_r = V_n / V_c;


    M = [           -lambda(1),                     0,                         0,       a(1),                       0,          0,                         0,             0,                              0;
         lambda(1) * V_c / V_n,            -lambda(1),                         0,          0,                       0,          0,                         0,             0,                              0;
                             0,                  beta, -(lambda(1) + lambda_ECF),          0,                       0,          0,                         0,             0,                     lambda_ECF;
                             0,                     0,     lambda(1) * V_n / V_c, -lambda(1),                       0,          0,                         0,             0,                              0;
                             0,                     0,                         0,          0,              -lambda(2),          0,                         0,          a(2),                              0;
                             0,                     0,                         0,          0,   lambda(2) * V_c / V_n, -lambda(2),                         0,             0,                              0;
                             0,                     0,                         0,          0,                       0,       beta, -(lambda(2) + lambda_ECF),             0,                     lambda_ECF;
                             0,                     0,                         0,          0,                       0,          0,     lambda(2) * V_n / V_c,    -lambda(2),                              0;
                             0,                     0,  lambda_ECF * V_c / V_ECF,          0,                       0,          0,  lambda_ECF * V_c / V_ECF,             0, -2 * lambda_ECF * V_c / V_ECF];


    z(:, j) = eig(M);
    if j == 2
        lambda_ECF = temp;
    end
end

disp(real(z(:, 1)));
[MY,IY] = max(real(z(:, 1)));
disp(IY);

disp(real(z(:, 1)));
[MN,IN] = max(real(z(:, 1)));
disp(IN);

% plot(real(z(:, 2)), imag(z(:, 2)), '*r');
% grid on;
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';

figure(1);

subplot(1,2,1)
plot(real(z(:, 1)), imag(z(:, 1)), '*r');
grid on;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title(['Eigenvalues of A with \lambda_{ECF} = ', num2str(lambda_ECF)])

subplot(1,2,2)
plot(real(z(:, 2)), imag(z(:, 2)), '*r');
grid on;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
title(['Eigenvalues of A with \lambda_{ECF} = ', num2str(0)])
%  
% solve for steady-state (a non-linear system)
function F = SolveSteady(x, n_cells, alpha, V_n, K, r, lambda, V_c, beta, lambda_ECF, V_ECF)

	F = zeros(size(x));
    
	for i = 1:n_cells
		F((i-1)*4 + 1) = ((alpha / V_n) * (K / (K + x((i-1)*4 + 4))).^r(i) - lambda(i) * x((i-1)*4 + 1));
		F((i-1)*4 + 2) = (lambda(i) * (V_n / V_c) * x((i-1)*4 + 1) - lambda(i) * x((i-1)*4 + 2));
		F((i-1)*4 + 3) = (beta * x((i-1)*4 + 2) - lambda(i) * x((i-1)*4 + 3) ...
						- lambda_ECF * (x((i-1)*4 + 3) - (V_ECF / V_c) * x(end)));
		F((i-1)*4 + 4) = (lambda(i) * (V_c / V_n) * x((i-1)*4 + 3) - lambda(i) * x((i-1)*4 + 4));
		F(end) = F(end) + lambda_ECF * ((V_c / V_ECF) * x((i-1)*4 + 3) - x(end));
    end
    
end



