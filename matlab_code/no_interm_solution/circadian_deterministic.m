clc;
clear

% the variables in the circadian oscillator are m_n, m_c, p_c and p_n,
% which are number of mRNA moleculs in nucleus, number of mRNA molecules in
% cytoplasm, number of this particular type of protein molecules in
% cytoplasm, and number of protein molecules in nuclues, respectively.

% the parameters that control the ocsillation are alpha, gamma_m, delta_m,
% beta, gamma_p, delta_p, which are 


% initialize variables

%{
% p_c starts from the highest
p_n = 1403;
p_c = 1670;
m_n = 33;
m_c = 48;
%}

p_n = 500;
p_c = 500;
m_n = 500;
m_c = 500;

%{
% p_c starts from the lowest
p_n = 925;
p_c = 730;
m_n = 32;
m_c = 21;
%}

lambda = 0.285599;
gamma_m = lambda;
delta_m = lambda;
delta_p = lambda;
gamma_p = lambda;
r = 5;
k = 200;
alpha = 1e5;
beta = 10;


dt = 0.01;
endtime = 300;
timevec = 0:dt:endtime;
plot_to_the_end = length(timevec);
% plot_to_the_end = int32(2 * length(timevec) / 3);


pn_mat = zeros(length(timevec), 1);
pc_mat = zeros(length(timevec), 1);
mn_mat = zeros(length(timevec), 1);
mc_mat = zeros(length(timevec), 1);

for t = 1:length(timevec)
    m_n_new = (alpha * (k/(k+p_n)).^r - gamma_m * m_n)*dt + m_n;
    m_c_new = (gamma_m * m_n - delta_m * m_c)*dt + m_c;
    p_c_new = (beta * m_c - gamma_p * p_c)*dt + p_c;
    p_n_new = (gamma_p * p_c - delta_p * p_n)*dt + p_n;
    
    m_n = m_n_new;
    m_c = m_c_new;
    p_c = p_c_new;
    p_n = p_n_new;
    
    
    pn_mat(t) = p_n;
    pc_mat(t) = p_c;
    mn_mat(t) = m_n;
    mc_mat(t) = m_c;
    
end


%{ 
%ONLY NEEDED TO ADJUST INIT VALUES
pn_mat(1925,:)
pc_mat(1925,:)
mn_mat(1925,:)
mc_mat(1925,:)
%}

figure(1);

subplot(2,2,1)
plot(timevec, pn_mat);
title('p_n')
grid;

%legend({'P_n'},'Location','southwest')

subplot(2,2,2)
plot(timevec, pc_mat);
%legend({'P_c'},'Location','southwest')
title('p_c')
grid;


subplot(2,2,3)
plot(timevec, mn_mat);
title('m_n')
grid;


subplot(2,2,4)
plot(timevec, mc_mat);
title('m_c')
grid;

%suptitle('A cellular circadian oscillator')


%legend({'M_n', 'M_c'},'Location','southwest')


