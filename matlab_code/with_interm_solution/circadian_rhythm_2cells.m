% the variables in the circadian oscillator are m_n, m_c, p_c and p_n,
% which are number of mRNA moleculs in nucleus, number of mRNA molecules in
% cytoplasm, number of this particular type of protein molecules in
% cytoplasm, and number of protein molecules in nuclues, respectively.

% the parameters that control the ocsillation are alpha, gamma_m, delta_m,
% beta, gamma_p, delta_p, which are 


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
K_ex = 200;
K_im = 60;
alpha = 0.1;
beta = 0.1;
alpha_ex = 0;
alpha_im = 0;
r_ex = 13;
r_im = 15;
V_n = 1e-5;
V_c = 5e-5;
%% initialize gamma_m = delta_m = delta_p = gamma_p = lambda as well as r
%  different for each cell
r = zeros(1, n_cells);
r(1) = 10;
r(2) = 10;
r(3:end) = 10;
lambda = zeros(1, n_cells);
lambda(1) = 3 * pi / 220;
lambda(2) = 2 * pi / 110;
lambda(3:end) = 2.5 * pi / 110;
%% set timevec and init arrays for holding data
dt = 0.01;
endtime = 2000;
timevec = 0:dt:endtime;
plot_to_the_end = length(timevec);
% plot_to_the_end = int32(2 * length(timevec) / 3);

pn_mat = zeros(plot_to_the_end, n_cells);
pc_mat = zeros(plot_to_the_end, n_cells);
mn_mat = zeros(plot_to_the_end, n_cells);
mc_mat = zeros(plot_to_the_end, n_cells);
pECF_mat = zeros(plot_to_the_end, 1);


%% main loop
for t = 1:length(timevec)
    p_ECF = all_values(end);
    d_pECF = 0;
    
    if t > int32(1 * length(timevec) / 4)
        alpha_ex = 2e4;
        alpha_im = 2e3;
    end
    
    for i = 1:n_cells
        m_n = all_values((i - 1) * 4 + 1);
        m_c = all_values((i - 1) * 4 + 2);
        p_c = all_values((i - 1) * 4 + 3);
        p_n = all_values((i - 1) * 4 + 4);
        
        m_n_new = ((alpha / V_n) * (K / (K + p_n)).^r(i) - lambda(i) * m_n)...
            * dt + m_n;
        m_c_new = (lambda(i) * (V_n / V_c) * m_n - lambda(i) * m_c) ...
            * dt + m_c;
        p_c_new = (beta * m_c - lambda(i) * p_c ... 
            - (alpha_ex / V_c) * (p_c / (K_ex + p_c)).^r_ex ...
            + (alpha_im / V_ECF) * (p_ECF / (K_im + p_ECF)).^r_im) ... 
            * dt + p_c;
        p_n_new = (lambda(i) * (V_c / V_n) * p_c - lambda(i) * p_n) ...
            * dt + p_n;
        d_pECF = (alpha_ex / V_c) * (p_c / (K_ex + p_c)).^r_ex ... 
            * dt + d_pECF;
        
        m_n_new = max(m_n_new, 0);
        m_c_new = max(m_c_new, 0);
        p_c_new = max(p_c_new, 0);
        p_n_new = max(p_n_new, 0);
        
        all_values((i - 1) * 4 + 1) = m_n_new;
        all_values((i - 1) * 4 + 2) = m_c_new;
        all_values((i - 1) * 4 + 3) = p_c_new;
        all_values((i - 1) * 4 + 4) = p_n_new;

        if t > length(timevec) - plot_to_the_end
            pn_mat(t - length(timevec) + plot_to_the_end, i) = p_n_new;
            pc_mat(t - length(timevec) + plot_to_the_end, i) = p_c_new;
            mn_mat(t - length(timevec) + plot_to_the_end, i) = m_n_new;
            mc_mat(t - length(timevec) + plot_to_the_end, i) = m_c_new;
        end
        
        d_pECF = - (alpha_im / V_ECF) * (p_ECF / (K_im + p_ECF)).^r_im ...
            * n_cells * dt + d_pECF;
        
        pECF_new = max(p_ECF + d_pECF, 0);
        all_values(end) = pECF_new;
        if t > length(timevec) - plot_to_the_end
            pECF_mat(t - length(timevec) + plot_to_the_end) = pECF_new;
        end
    end
end

%% Plot
% plot each cell
for j = 1:n_cells
    figure(j);
    subplot(2,2,1)
    plot(pn_mat(:, j));
    title('p_n')
    grid;


    subplot(2,2,2)
    plot(pc_mat(:, j));
    title('p_c')
    grid;


    subplot(2,2,3)
    plot(mn_mat(:, j));
    title('m_n')
    grid;


    subplot(2,2,4)
    plot(mc_mat(:, j));
    title('m_c')
    grid;
end

% plot ECF
figure(j + 1);
plot(pECF_mat);
legend({'P_{ECF}'},'Location','southwest')
grid;