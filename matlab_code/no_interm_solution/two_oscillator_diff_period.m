clear; clc; close all;

% USE THE SIMPLEST WAY TO CONNECT TWO CELLS BY DIRECTLY ALLOW PROTEINS TO
% BE EXCHANGED.

% WE ASSUME THE EQUATION THAT GOVERNING THIS IS 
% dp1/dt = -C(p1-p2)
% dp2/dt = -dp1/dt


% the variables in the circadian oscillator are m_n, m_c, p_c and p_n,
% which are number of mRNA moleculs in nucleus, number of mRNA molecules in
% cytoplasm, number of this particular type of protein molecules in
% cytoplasm, and number of protein molecules in nuclues, respectively.
% the parameters that control the ocsillation are alpha, gamma_m, delta_m,
% beta, gamma_p, delta_p, which are 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variables and parameters for CELL 1
p1_n = 925;
p1_c = 730;
m1_n = 32;
m1_c = 21;
lambda1 = 2*pi/22;
gamma1_m = lambda1;
delta1_m = lambda1;
delta1_p = lambda1;
gamma1_p = lambda1;
r1 = 5;
k1 = 200;
alpha1 = 1e5;
beta1 = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variables and parameters for CELL 2
p2_n = 1403;
p2_c = 1670;
m2_n = 33;
m2_c = 48;

lambda2 = 2*pi/22*0.8;
gamma2_m = lambda2;
delta2_m = lambda2;
delta2_p = lambda2;
gamma2_p = lambda2;
r2 = 5;
k2 = 200;
alpha2 = 1e5;
beta2 = 10;

% Diffusion Constant don't know how large it is
F = .01;

% 0.01, 0.04, 0.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMPLEST CASE: THE TWO CELLS HAVE SAME PERIOD BUT JUST DIFFERENT INITIAL
% STATES


dt = 0.01;
endtime = 2000;
timevec = 0:dt:endtime;
plot_to_the_end = length(timevec);


p1n_mat = zeros(length(timevec), 1);
p1c_mat = zeros(length(timevec), 1);
m1n_mat = zeros(length(timevec), 1);
m1c_mat = zeros(length(timevec), 1);

p2n_mat = zeros(length(timevec), 1);
p2c_mat = zeros(length(timevec), 1);
m2n_mat = zeros(length(timevec), 1);
m2c_mat = zeros(length(timevec), 1);

diffusion_mat = zeros(length(timevec), 1);

for t = 1:length(timevec)
    
    % reactions in cell 1
    m1_n_new = (alpha1 * (k1/(k1+p1_n)).^r1 - gamma1_m * m1_n)*dt + m1_n;
    m1_c_new = (gamma1_m * m1_n - delta1_m * m1_c)*dt + m1_c;
    p1_c_new = (beta1 * m1_c - gamma1_p * p1_c)*dt + p1_c;
    p1_n_new = (gamma1_p * p1_c - delta1_p * p1_n)*dt + p1_n;
    
    m1_n = m1_n_new;
    m1_c = m1_c_new;
    p1_c = p1_c_new;
    p1_n = p1_n_new;
    
    
    % reactions in cell 2
    m2_n_new = (alpha2 * (k2/(k2+p2_n)).^r2 - gamma2_m * m2_n)*dt + m2_n;
    m2_c_new = (gamma2_m * m2_n - delta2_m * m2_c)*dt + m2_c;
    p2_c_new = (beta2 * m2_c - gamma2_p * p2_c)*dt + p2_c;
    p2_n_new = (gamma2_p * p2_c - delta2_p * p2_n)*dt + p2_n;
    
    m2_n = m2_n_new;
    m2_c = m2_c_new;
    p2_c = p2_c_new;
    p2_n = p2_n_new;
    
    % interactions between cells
    % this is given by 
    % d{p1_c}/dt = - F * (p1_c - p2_c)
    % d{p2_c}/dt = - d{p1_c}/dt
    p1_c_new = p1_c - F * (p1_c - p2_c) * dt;
    p2_c_new = p2_c + F * (p1_c - p2_c) * dt;
    
    p1_c = p1_c_new;
    p2_c = p2_c_new;
    
    p1n_mat(t) = p1_n;
    p1c_mat(t) = p1_c;
    m1n_mat(t) = m1_n;
    m1c_mat(t) = m1_c;
    
    p2n_mat(t) = p2_n;
    p2c_mat(t) = p2_c;
    m2n_mat(t) = m2_n;
    m2c_mat(t) = m2_c;
    
    diffusion_mat(t) = F * (p1_c - p2_c) * dt;
    
end


%{ 
%ONLY NEEDED TO ADJUST INIT VALUES
pn_mat(1925,:)
pc_mat(1925,:)
mn_mat(1925,:)
mc_mat(1925,:)
%}



figure(1);

subplot(3,1,1)
plot(timevec, p1c_mat);
hold on;
grid;
plot(timevec, p2c_mat);
hold off;
title(['p^1_c and p^2_c.',' At F = ', num2str(F)])
legend('p^1_c', 'p^2_c','Location','northeast')

subplot(3,1,2)
plot(timevec, m1n_mat)
hold on;
plot(timevec, m2n_mat)
legend('m^1_n', 'm^2_n','Location','northeast')
hold off;
grid;
title(['m^1_n and m^2_n.',' At F = ', num2str(F)])

subplot(3,1,3)
plot(timevec, diffusion_mat)
title(['Diffusion of p_c at F = ', num2str(F)])
grid;
title(['Flux of p_c between the two cells. ',' At F = ', num2str(F)])
set(gcf,'Position',[100 100 1000 500])


%{
figure(1);

subplot(3,2,1)
plot(timevec, p1c_mat);
title('p^1_c')
grid;

subplot(3,2,2)
plot(timevec, p2c_mat);
title('p^2_c')
grid;

subplot(3,2,3)
plot(timevec, m1n_mat);
title('m^1_n')
grid;

subplot(3,2,4)
plot(timevec, m2n_mat);
title('m^2_n')
grid;

subplot(3,2,[5,6])
plot(timevec, diffusion_mat)
title(['Diffusion of p_c at F = ', num2str(F)])
grid;
title(['Flux of p_c between the two cells. ','F = ', num2str(F)])

%suptitle(['F = ', num2str(F)])

figure(2)
subplot(2,1,1)
plot(timevec, m1n_mat)
hold;
plot(timevec, m2n_mat)
legend('m^1_n', 'm^2_n','Location','northeast')
%title(
subplot(2,1,2)
plot(timevec, diffusion_mat)
legend(['Flux of p_c between the two cells. ','F = ', num2str(F)])

%}
%{
subplot(2,2,2)
plot(timevec, p1c_mat);
title('p1_c')
grid;



subplot(2,2,3)
plot(timevec, m1n_mat);
title('m1_n')
grid;


subplot(2,2,4)
plot(timevec, m1c_mat);
title('m1_c')
grid;
%}


