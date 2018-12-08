clear;
clc;
% the code below illustrates discrete/stochastic version of circadian
% rhythm of a cell


% parameters 

r = 5;
k = 200;
alpha = 1e5;
beta = 10;
lambda = 0.285599;
gamma_m = lambda;
gamma_p = lambda;
delta_m = lambda;
delta_p = lambda;


% initial values
p_n = 925;
p_c = 730;
m_n = 32;
m_c = 21;


t = 0;
tmax = 40;

% for plotting only
count_pn = 1;
count_pc = 1;
count_mn = 1;
count_mc = 1;
plot_time_vec_pn(count_pn) = 0;
plot_time_vec_pc(count_pc) = 0;
plot_time_vec_mn(count_mn) = 0;
plot_time_vec_mc(count_mc) = 0;

pn_mat(count_pn) = p_n;
pc_mat(count_pc) = p_c;
mn_mat(count_mn) = m_n;
mc_mat(count_mc) = m_c;

pn_change = 0;
pc_change = 0;
mn_change = 0;
mc_change = 0;


while t < tmax
    rate = rates(m_n, m_c, p_c, p_n);

    % rate is a 6*1 matrix that denotes the possibility of each reaction
    % happening
    
    T = - log(rand(1,6))./rate;
    % T is a 6*1 matrix that denotes the time it takes for each of the
    % reaction to happen. We take the idea that only the one with the
    % smallest time indeed happens at this stage.
    
    [Tmin, kmin] = min(T);
    % Tmin is the time it takes for the "happened" reaction to happen. kmin
    % is the position (1,2,3,4,5, or 6) indicating which reacation happens
    
    switch kmin
        case 1
            % transcription happens
            m_n = m_n + 1;
            
            % for plotting only
            mn_change = 1;
        case 2
            % mRNA enters cytoplasm from necleus
            m_n = m_n - 1;
            m_c = m_c + 1;
            
            % for plotting only
            mn_change = 1;
            mc_change = 1;
        case 3
            % translation happens
            p_c = p_c + 1;
            
            % for plotting only
            pc_change = 1;
        case 4 
            % protein enters nucleus from cytoplasm
            p_c = p_c - 1;
            p_n = p_n + 1;
            
            % for plotting only
            pc_change = 1;
            pn_change = 1;
        case 5
            % degradation of mRNA
            m_c = m_c - 1;
            
            % for plotting only
            mc_change = 1;
        case 6
            % degradation of protein
            p_n = p_n - 1;
            
            % for plotting only
            pn_change = 1;
    end
    
 
    t = t + Tmin;
    
    if pc_change == 1
        count_pc = count_pc + 1;
        plot_time_vec_pc(count_pc) = t;
        pc_mat(count_pc) = p_c;
        pc_change = 0;
    end
    
    if pn_change == 1
        count_pn = count_pn + 1;
        plot_time_vec_pn(count_pn) = t;
        pn_mat(count_pn) = p_n;
        pn_change = 0;
    end
    
    if mc_change == 1
        count_mc = count_mc + 1;
        plot_time_vec_mc(count_mc) = t;
        mc_mat(count_mc) = m_c;
        mc_change = 0;
    end
    
    if mn_change == 1
        count_mn = count_mn + 1;
        plot_time_vec_mn(count_mn) = t;
        mn_mat(count_mn) = m_n;
        mn_change = 0;
    end
    

end

figure(1)

% neglect the last point for purpose of drawing
subplot(2,2,1)
stairs(plot_time_vec_pn(1:length(pn_mat)-1), pn_mat(1:length(pn_mat)-1))
title('p_n')


subplot(2,2,2)
stairs(plot_time_vec_pc(1:length(pc_mat)-1), pc_mat(1:length(pc_mat)-1))
title('p_c')


subplot(2,2,3)
stairs(plot_time_vec_mn(1:length(mn_mat)-1), mn_mat(1:length(mn_mat)-1))
title('m_n')


subplot(2,2,4)
stairs(plot_time_vec_mc(1:length(mc_mat)-1), mc_mat(1:length(mc_mat)-1))
title('m_c')


function rate= rates(m_n, m_c, p_c, p_n)
    %[alpha, beta, k, r, gamma_m, gamma_p, delta_m, delta_p] = params;
    %alpha
    %rate= zeros(6);
    r = 5;
    k = 200;
    alpha = 1e5;
    beta = 10;
    lambda = 2*pi / 22;
    gamma_m = lambda;
    gamma_p = lambda;
    delta_m = lambda;
    delta_p = lambda;
    
    rate(1) = alpha * (k/(k+p_n)).^r;
    rate(2) = gamma_m * m_n;
    rate(3) = beta * m_c;
    rate(4) = gamma_p * p_c;
    rate(5) = delta_m * m_c;
    rate(6) = delta_p * p_n;
end