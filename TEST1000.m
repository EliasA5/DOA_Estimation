close all
clear all
clc
%--------------------------------------------------------------------------
K_1 = 15;
K_3 = 1;
alpha = pi/4;
v_0 = 9600;  % P waves velocities can be from 6 Km/sec to 11 Km/sec
SNR = 10;
sigma_source = 10;
sigma_noise = 1;  % for white and colored noise
data = load("data.mat");
r_k = data.r_m;   % r_m is  loaded form the data matrix we extrcted from the getData python script
M = 100;   % for these simlutaions we decided to use one segement
P = 1;
f = 3;      % the omega_m is constant at this point for all the frquencies
%--------------------------------------------------------------------------
Tests = 10;

theta_og = -pi +2*pi/Tests : 2*pi/Tests : pi;
theta_0 = pi/4; % starting estimate at 45 deg

iters = 400;

theta_white_f = zeros(Tests,1);
theta_colored_f = zeros(Tests,1);


for t = 1 : Tests

    step_size_white = 1;
    step_size_colored = 1;
    gamma = 0.85;
    
    theta_colored = zeros(1,iters);
    theta_white = zeros(1,iters);
    
    theta_white(1) = theta_0;
    theta_colored(1) = theta_0;
    
    [X_colored,s_colored,Rv_colored,~] = synData(r_k, theta_og(t), alpha, v_0, sigma_source, sigma_noise, M, 'colored', f, K_1, K_3, P);
    [X_white,s_white,Rv_white,~] = synData(r_k, theta_og(t), alpha, v_0, sigma_source, sigma_noise, M, 'white', f, K_1, K_3, P);
    
    a_colored = zeros(3*K_3 + K_1 , 1);
    a_white = zeros(3*K_3 + K_1 , 1);
    
    da_colored = zeros(3*K_3 + K_1 , 1);
    da_white = zeros(3*K_3 + K_1 , 1);
    
    Rv_white_inv = pinv(Rv_white);
    Rv_colored_inv = pinv(Rv_colored);
    
    X_w_white = squeeze(sum(X_white,1));
    X_w_colored = squeeze(sum(X_colored,1));
    
    sin_alpha = sin(alpha);
    cos_alpha = cos(alpha);

    pdf_white = zeros(1,iters);
    pdf_colored = zeros(1,iters);

    for i = 1 : iters - 1
    
        FIM_white = 0;
        FIM_colored = 0;
    
        df_theta_white = 0;
        df_theta_colored = 0;
    
        % calculate the a and a derivative for the colored and white
        sin_theta_white = sin(theta_white(i));
        sin_theta_colored = sin(theta_colored(i));
    
        cos_theta_white = cos(theta_white(i));
        cos_theta_colored = cos(theta_colored(i));
    
        e_u_white = [sin_theta_white * sin_alpha ; cos_theta_white * sin_alpha ; cos_alpha];
        e_u_colored = [sin_theta_colored * sin_alpha ; cos_theta_colored * sin_alpha ; cos_alpha];
    
        for j = 1 : 3 : 3 * K_3
    
            r_k_i = r_k(j,:);
    
            tau_3D_white = r_k_i * e_u_white/v_0;
            tau_3D_colored = r_k_i * e_u_colored/v_0;
    
            tau_3D_der_white = (1/v_0) * (r_k_i(1) * cos_theta_white * sin_alpha - r_k_i(2) * sin_theta_white * sin_alpha);
            tau_3D_der_colored = (1/v_0) * (r_k_i(1) * cos_theta_colored * sin_alpha - r_k_i(2) * sin_theta_colored * sin_alpha);
    
            e_i_white = exp(-1i * f * tau_3D_white);
            e_i_colored = exp(-1i * f * tau_3D_colored);
    
            da_white(j) = sin_alpha * e_i_white * (cos_theta_white - 1i * f * sin_theta_white * tau_3D_der_white);
            da_white(j+1) = sin_alpha * e_i_white * (-sin_theta_white - 1i * f * cos_theta_white * tau_3D_der_white);
            da_white(j+2) = cos_alpha * e_i_white * (-1i * f * tau_3D_der_white);
        
            da_colored(j) = sin_alpha * e_i_colored * (cos_theta_colored - 1i * f * sin_theta_colored * tau_3D_der_colored);
            da_colored(j+1) = sin_alpha * e_i_colored * (-sin_theta_colored - 1i * f * cos_theta_colored * tau_3D_der_colored);
            da_colored(j+2) = cos_alpha * e_i_colored * (-1i * f * tau_3D_der_colored);
    
            a_white(j) = sin_theta_white * sin_alpha * e_i_white;
            a_white(j+1) = cos_theta_white * sin_alpha * e_i_white;
            a_white(j+2) = cos_alpha * e_i_white;
        
            a_colored(j) = sin_theta_colored * sin_alpha * e_i_colored;
            a_colored(j+1) = cos_theta_colored * sin_alpha * e_i_colored;
            a_colored(j+2) = cos_alpha * e_i_colored;
    
        end
    
        for j = 3 * K_3 + 1 : 3 * K_3 + K_1
    
            r_k_i = r_k(j,:);
    
            tau_1D_white = r_k_i * e_u_white/v_0;
            tau_1D_colored = r_k_i * e_u_colored/v_0;
    
            tau_1D_der_white = (1/v_0) * (r_k_i(1) * cos_theta_white * sin_alpha - r_k_i(2) * sin_theta_white * sin_alpha);
            tau_1D_der_colored = (1/v_0) * (r_k_i(1) * cos_theta_colored * sin_alpha - r_k_i(2) * sin_theta_colored * sin_alpha);
    
            e_i_white = exp(-1i * f * tau_1D_white);
            e_i_colored = exp(-1i * f * tau_1D_colored);
    
            da_white(j) = cos_alpha * e_i_white * (-1i * f * tau_1D_der_white);
            da_colored(j) = cos_alpha * e_i_colored * (-1i * f * tau_1D_der_colored);
    
            a_white(j) = cos_alpha * e_i_white;
            a_colored(j) = cos_alpha * e_i_colored;
    
        end
    
        for j = 1 : M
    
            mue_white = a_white * s_white(j);
            mue_colored = a_colored * s_colored(j);
    
            dmue_white = da_white * s_white(j);
            dmue_colored = da_colored * s_colored(j);
    
            x_i_white = X_w_white(:,j);
            x_i_colored = X_w_colored(:,j);
    
            df_theta_white = df_theta_white + (x_i_white - mue_white)' * ...
            Rv_white_inv * dmue_white + (dmue_white' * Rv_white_inv * (x_i_white - mue_white)).';
    
            df_theta_colored = df_theta_colored + (x_i_colored - mue_colored)' * ...
            Rv_colored_inv * dmue_colored + (dmue_colored' * Rv_colored_inv * (x_i_colored - mue_colored)).';
    
            FIM_white = FIM_white + real(dmue_white' * Rv_white_inv * dmue_white);
            FIM_colored = FIM_colored + real(dmue_colored' * Rv_colored_inv * dmue_colored);
    
        end
    
        FIM_white = FIM_white * 2;
        FIM_colored = FIM_colored * 2;

        pdf_white(i+1) = log_likelihood(X_w_white,mue_white,Rv_white_inv,M);
        pdf_colored(i+1) = log_likelihood(X_w_colored,mue_colored,Rv_colored_inv,M);
    
        theta_white(i+1) = wrapToPi(real(theta_white(i) + step_size_white * FIM_white ^ (-1) * df_theta_white));
        theta_colored(i+1) = wrapToPi(real(theta_colored(i) + step_size_colored * FIM_colored ^ (-1) * df_theta_colored));

        if (pdf_white(i+1) < pdf_white(i))
            step_size_white = step_size_white * gamma;
        end

        if (pdf_colored(i+1) < pdf_colored(i))
            step_size_colored = step_size_colored * gamma;
        end
    end
    
    theta_white_f(t) = (theta_white(end));
    theta_colored_f(t) = (theta_colored(end));

%     figure;
%     hold on
%     plot(theta_colored)
%     plot(theta_white)
%     hold off
%     legend('colored' , white)

end

bias_white = wrapToPi((theta_white_f - theta_og.'));
bias_colored = wrapToPi((theta_colored_f - theta_og.'));

% MSP_white = wrapToPi((theta_white_f - theta_og.') .^ 2);
% MSP_colored = wrapToPi((theta_colored_f - theta_og.') .^ 2);

MSP_white = (cos(theta_og.') - cos(theta_white_f) .^ 2) + (sin(theta_og.') - sin(theta_white_f) .^ 2);
MSP_colored = (cos(theta_og.') - cos(theta_colored_f) .^ 2) + (sin(theta_og.') - sin(theta_colored_f) .^ 2);


figure;
hold on;
plot(theta_og,'*')
plot(theta_white_f,'o')
stem(theta_colored_f)
legend('Original', 'White Estimate', 'Colored Estimate')
xlabel('Test')
ylabel('\theta')
hold off

figure;
hold on;
stem(bias_white)
stem(bias_colored)
legend('White Estimate','Colored Estimate')
xlabel('Test')
ylabel('Bias')
hold off

figure;
hold on;
stem(MSP_white)
stem(MSP_colored)
legend('White Estimate','Colored Estimate')
xlabel('Test')
ylabel('MSE')
hold off

%%
angles = linspace(0, 2*pi, 500);
radius = 20;
CenterX = 50;
CenterY = 40;
x = radius * cos(angles) + CenterX;
y = radius * sin(angles) + CenterY;

x_og = radius * cos(theta_og) + CenterX;
x_white = radius * cos(theta_white_f) + CenterX;
x_colored = radius * cos(theta_colored_f) + CenterX;

y_og = radius * sin(theta_og) + CenterY;
y_white = radius * sin(theta_white_f) + CenterY;
y_colored = radius * sin(theta_colored_f) + CenterY;

figure
hold on
grid on
axis equal;
plot(x, y, 'k-', 'LineWidth', 1);
plot(x_og, y_og, 'rx','LineWidth', 2)
plot(x_white, y_white, 'bo','LineWidth', 1.5)
plot(x_colored, y_colored, 'go','LineWidth', 1.5)
xlabel('X', 'FontSize', 10);
ylabel('Y', 'FontSize', 10);
legend('','Original','White','Colored')



