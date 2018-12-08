clear;
clc;
% the code below illustrates discrete/stochastic version of circadian
% rhythm of two cells


% cell 1
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

% cell 2
p2_n = 1403;
p2_c = 1670;
m2_n = 33;
m2_c = 48;

lambda2 = 2*pi/22;
gamma2_m = lambda2;
delta2_m = lambda2;
delta2_p = lambda2;
gamma2_p = lambda2;
r2 = 5;
k2 = 200;
alpha2 = 1e5;
beta2 = 10;

% Diffusion Constant
F = 0.08;


t = 0;
tmax = 1000;

% for plotting only
count_p1n = 1;
count_p1c = 1;
count_m1n = 1;
count_m1c = 1;

count_p2n = 1;
count_p2c = 1;
count_m2n = 1;
count_m2c = 1;

count_diffusion = 1;

plot_time_vec_p1n(count_p1n) = 0;
plot_time_vec_p1c(count_p1c) = 0;
plot_time_vec_m1n(count_m1n) = 0;
plot_time_vec_m1c(count_m1c) = 0;

plot_time_vec_p2n(count_p2n) = 0;
plot_time_vec_p2c(count_p2c) = 0;
plot_time_vec_m2n(count_m2n) = 0;
plot_time_vec_m2c(count_m2c) = 0;

plot_time_vec_diffusion(count_diffusion) = 0;

p1n_mat(count_p1n) = p1_n;
p1c_mat(count_p1c) = p1_c;
m1n_mat(count_m1n) = m1_n;
m1c_mat(count_m1c) = m1_c;

p2n_mat(count_p2n) = p2_n;
p2c_mat(count_p2c) = p2_c;
m2n_mat(count_m2n) = m2_n;
m2c_mat(count_m2c) = m2_c;

diffusion_mat(count_diffusion) = 0;

p1n_change = 0;
p1c_change = 0;
m1n_change = 0;
m1c_change = 0;

p2n_change = 0;
p2c_change = 0;
m2n_change = 0;
m2c_change = 0;

is_diffusion = 0;
flag1 = 0;
flag2 = 0;


while t < tmax
    rate = rates(m1_n, m1_c, p1_c, p1_n, m2_n, m2_c, p2_c, p2_n, F);

    % rate is a 6*1 matrix that denotes the possibility of each reaction
    % happening
    
    T = - log(rand(1,13))./rate;
    % T is a 6*1 matrix that denotes the time it takes for each of the
    % reaction to happen. We take the idea that only the one with the
    % smallest time indeed happens at this stage.
    
    [Tmin, kmin] = min(T);
    % Tmin is the time it takes for the "happened" reaction to happen. kmin
    % is the position (1,2,3,4,5, or 6) indicating which reacation happens
    
    switch kmin
        case 1
            % transcription in cell 1 happens
            m1_n = m1_n + 1;
            
            % for plotting only
            m1n_change = 1;
        case 2
            % mRNA enters cytoplasm from necleus in cell 1
            m1_n = m1_n - 1;
            m1_c = m1_c + 1;
            
            % for plotting only
            m1n_change = 1;
            m1c_change = 1;
        case 3
            % translation happens in cell 1
            p1_c = p1_c + 1;
            
            % for plotting only
            p1c_change = 1;
        case 4 
            % protein enters nucleus from cytoplasm in cell 1
            p1_c = p1_c - 1;
            p1_n = p1_n + 1;
            
            % for plotting only
            p1c_change = 1;
            p1n_change = 1;
        case 5
            % degradation of mRNA in cell 1
            m1_c = m1_c - 1;
            
            % for plotting only
            m1c_change = 1;
        case 6
            % degradation of protein in cell 1
            p1_n = p1_n - 1;
            
            % for plotting only
            p1n_change = 1;
            
        case 7
            % transcription in cell 2 happens
            m2_n = m2_n + 1;
            
            % for plotting only
            m2n_change = 1;
        case 8
            % mRNA enters cytoplasm from necleus in cell 2
            m2_n = m2_n - 1;
            m2_c = m2_c + 1;
            
            % for plotting only
            m2n_change = 1;
            m2c_change = 1;
        case 9
            % translation happens in cell 2
            p2_c = p2_c + 1;
            
            % for plotting only
            p2c_change = 1;
        case 10
            % protein enters nucleus from cytoplasm in cell 2
            p2_c = p2_c - 1;
            p2_n = p2_n + 1;
            
            % for plotting only
            p2c_change = 1;
            p2n_change = 1;
        case 11
            % degradation of mRNA in cell 2
            m2_c = m2_c - 1;
            
            % for plotting only
            m2c_change = 1;
        case 12
            % degradation of protein in cell 2
            p2_n = p2_n - 1;
            
            % for plotting only
            p2n_change = 1;
            
        case 13
            % diffusion between cells
            if p1_c > p2_c
                % should flow from cell 1 to cell 2
                p1_c = p1_c - 1;
                p2_c = p2_c + 1;
                
                % a flag to denote which direction the diffusion is
                flag1 = 1;
            else
                % should flow from cell 2 to cell 1
                p2_c = p2_c - 1;
                p1_c = p1_c + 1;
                
                % a flag to denote which direction the diffusion is
                flag2 = 1;
            end
            
            % for plotting only
            p1c_change = 1;
            p2c_change = 1;
            is_diffusion = 1;
    end
    
 
    t = t + Tmin;
    
    if p1c_change == 1
        count_p1c = count_p1c + 1;
        plot_time_vec_p1c(count_p1c) = t;
        p1c_mat(count_p1c) = p1_c;
        p1c_change = 0;
    end
    
    if p1n_change == 1
        count_p1n = count_p1n + 1;
        plot_time_vec_p1n(count_p1n) = t;
        p1n_mat(count_p1n) = p1_n;
        p1n_change = 0;
    end
    
    if m1c_change == 1
        count_m1c = count_m1c + 1;
        plot_time_vec_m1c(count_m1c) = t;
        m1c_mat(count_m1c) = m1_c;
        m1c_change = 0;
    end
    
    if m1n_change == 1
        count_m1n = count_m1n + 1;
        plot_time_vec_m1n(count_m1n) = t;
        m1n_mat(count_m1n) = m1_n;
        m1n_change = 0;
    end
    
    
    if p2c_change == 1
        count_p2c = count_p2c + 1;
        plot_time_vec_p2c(count_p2c) = t;
        p2c_mat(count_p2c) = p2_c;
        p2c_change = 0;
    end
    
    if p2n_change == 1
        count_p2n = count_p2n + 1;
        plot_time_vec_p2n(count_p2n) = t;
        p2n_mat(count_p2n) = p2_n;
        p2n_change = 0;
    end
    
    if m2c_change == 1
        count_m2c = count_m2c + 1;
        plot_time_vec_m2c(count_m2c) = t;
        m2c_mat(count_m2c) = m2_c;
        m2c_change = 0;
    end
    
    if m2n_change == 1
        count_m2n = count_m2n + 1;
        plot_time_vec_m2n(count_m2n) = t;
        m2n_mat(count_m2n) = m2_n;
        m2n_change = 0;
    end
    
    if is_diffusion == 1
        count_diffusion = count_diffusion + 1;
        plot_time_vec_diffusion(count_diffusion) = t;
        if flag1 == 1
            diffusion_mat(count_diffusion) = 1;
        elseif flag2 == 1
            diffusion_mat(count_diffusion) = -1;
        end
        flag1 = 0; flag2 = 0; is_diffusion = 0;
    end

end

%{
figure(1)

% neglect the last point for purpose of drawing
subplot(2,2,1)
stairs(plot_time_vec_p1n(1:length(p1n_mat)-1), p1n_mat(1:length(p1n_mat)-1))
title('p^1_n')


subplot(2,2,2)
stairs(plot_time_vec_p1c(1:length(p1c_mat)-1), p1c_mat(1:length(p1c_mat)-1))
title('p^1_c')


subplot(2,2,3)
stairs(plot_time_vec_m1n(1:length(m1n_mat)-1), m1n_mat(1:length(m1n_mat)-1))
title('m^1_n')


subplot(2,2,4)
stairs(plot_time_vec_m1c(1:length(m1c_mat)-1), m1c_mat(1:length(m1c_mat)-1))
title('m^1_c')

figure(2)

% neglect the last point for purpose of drawing
subplot(2,2,1)
stairs(plot_time_vec_p2n(1:length(p2n_mat)-1), p2n_mat(1:length(p2n_mat)-1))
title('p^2_n')


subplot(2,2,2)
stairs(plot_time_vec_p2c(1:length(p2c_mat)-1), p2c_mat(1:length(p2c_mat)-1))
title('p^2_c')


subplot(2,2,3)
stairs(plot_time_vec_m2n(1:length(m2n_mat)-1), m2n_mat(1:length(m2n_mat)-1))
title('m^2_n')


subplot(2,2,4)
stairs(plot_time_vec_m2c(1:length(m2c_mat)-1), m2c_mat(1:length(m2c_mat)-1))
title('m^2_c')
%}


figure(3);

subplot(3,1,1)
stairs(plot_time_vec_p1c(1:length(p1c_mat)-1), p1c_mat(1:length(p1c_mat)-1));
hold on;
grid;
stairs(plot_time_vec_p2c(1:length(p2c_mat)-1), p2c_mat(1:length(p2c_mat)-1));
hold off;
title(['p^1_c and p^2_c.',' At F = ', num2str(F)])
legend('p^1_c', 'p^2_c','Location','northeast')

subplot(3,1,2)
stairs(plot_time_vec_m1n(1:length(m1n_mat)-1), m1n_mat(1:length(m1n_mat)-1));
hold on;
stairs(plot_time_vec_m2n(1:length(m2n_mat)-1), m2n_mat(1:length(m2n_mat)-1));
legend('m^1_n', 'm^2_n','Location','northeast')
hold off;
grid;
title(['m^1_n and m^2_n.',' At F = ', num2str(F)])

subplot(3,1,3)
plot(plot_time_vec_diffusion(1:length(diffusion_mat-1)), diffusion_mat(1:length(diffusion_mat-1)), 'o');
ylim([-2 2]);
title(['Diffusion of p_c at F = ', num2str(F)])
grid;
title(['Flux of p_c between the two cells. ',' At F = ', num2str(F)])

set(gcf,'Position',[100 100 1000 500])



function rate = rates(m1_n, m1_c, p1_c, p1_n, m2_n, m2_c, p2_c, p2_n, F)
    
    % cell 1
    lambda1 = 2*pi/22;
    gamma1_m = lambda1;
    delta1_m = lambda1;
    delta1_p = lambda1;
    gamma1_p = lambda1;
    r1 = 5;
    k1 = 200;
    alpha1 = 1e5;
    beta1 = 10;

    % cell 2

    lambda2 = 2*pi/22;
    gamma2_m = lambda2;
    delta2_m = lambda2;
    delta2_p = lambda2;
    gamma2_p = lambda2;
    r2 = 5;
    k2 = 200;
    alpha2 = 1e5;
    beta2 = 10;

    
    
    rate(1) = alpha1 * (k1/(k1+p1_n)).^r1;
    rate(2) = gamma1_m * m1_n;
    rate(3) = beta1 * m1_c;
    rate(4) = gamma1_p * p1_c;
    rate(5) = delta1_m * m1_c;
    rate(6) = delta1_p * p1_n;
    
    rate(7) = alpha2 * (k2/(k2+p2_n)).^r2;
    rate(8) = gamma2_m * m2_n;
    rate(9) = beta2 * m2_c;
    rate(10) = gamma2_p * p2_c;
    rate(11) = delta2_m * m2_c;
    rate(12) = delta2_p * p2_n;
    
    rate(13) = F * abs(p1_c - p2_c);
    
end