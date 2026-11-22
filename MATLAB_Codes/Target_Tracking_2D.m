% Scenario 1: 2D Kalman Filter Tracking (No GPS Dropouts) - Tuned + LS Init
clear; clc; close all;

% Settings
dt = 0.1;
N  = 300;
t  = (0:N-1)*dt;

x0_true = [0; 0; 1.2; 0.7];   % [x;y;vx;vy]

% --- TRUE system noise (used to generate ground truth) ---
sigma_a_true = 0.5;   % (m/s^2) true acceleration disturbance

% --- GPS noise (used to generate measurements) ---
sigma_z_true = 0.5;   % (m) true GPS std

% --- FILTER tuning knobs ---
sigma_a_filt = 0.5;   % (m/s^2) process noise assumed by KF
sigma_z_filt = 0.5;   % (m) measurement noise assumed by KF

rng(7); % repeatable

% Model (Constant Velocity)
F = [1 0 dt 0;
     0 1 0 dt;
     0 0 1  0;
     0 0 0  1];

H = [1 0 0 0;
     0 1 0 0];

% Q uses FILTER sigma_a_filt (NOT sigma_a_true)
Q = sigma_a_filt^2 * [dt^4/4   0      dt^3/2  0;
                      0      dt^4/4   0      dt^3/2;
                      dt^3/2  0      dt^2    0;
                      0      dt^3/2   0      dt^2];

% R matches FILTER sigma_z_filt
R = sigma_z_filt^2 * eye(2);

% True trajectory (uses TRUE sigma_a_true)
Xtrue = zeros(4, N);
Xtrue(:,1) = x0_true;

for k = 2:N
    ax = sigma_a_true * randn();
    ay = sigma_a_true * randn();

    Xtrue(1,k) = Xtrue(1,k-1) + Xtrue(3,k-1)*dt + 0.5*ax*dt^2;
    Xtrue(2,k) = Xtrue(2,k-1) + Xtrue(4,k-1)*dt + 0.5*ay*dt^2;
    Xtrue(3,k) = Xtrue(3,k-1) + ax*dt;
    Xtrue(4,k) = Xtrue(4,k-1) + ay*dt;
end

% GPS measurements (always available, uses TRUE sigma_z_true)
Z = H*Xtrue + sigma_z_true*randn(2, N);

% Kalman Filter
xhat = zeros(4, N);
P    = zeros(4,4,N);
I = eye(4);

% Improved initialization using Least Squares on first M points
% This reduces the "wiggle" caused by estimating v from only 2 noisy points

M = 10;                         % first 10 samples (1 second since dt=0.1)
tt = (0:M-1)' * dt;
A  = [ones(M,1) tt];            % pos(t) = p0 + v*t

bx = A \ Z(1,1:M)';             % [x0_est; vx_est]
by = A \ Z(2,1:M)';             % [y0_est; vy_est]

xhat(:,1) = [bx(1); by(1); bx(2); by(2)];

% Initial covariance:
% position ~ GPS noise, velocity uncertainty larger 
P(:,:,1)  = diag([sigma_z_filt^2, sigma_z_filt^2, 20^2, 20^2]);

for k = 2:N
    % Predict
    x_pred = F*xhat(:,k-1);
    P_pred = F*P(:,:,k-1)*F' + Q;

    % Update
    y = Z(:,k) - H*x_pred;
    S = H*P_pred*H' + R;
    K = P_pred*H' / S;

    xhat(:,k) = x_pred + K*y;
    P(:,:,k)  = (I - K*H) * P_pred;
end

% Evaluation (RMSE)
err_meas = vecnorm(Xtrue(1:2,:) - Z, 2, 1);
err_kf   = vecnorm(Xtrue(1:2,:) - xhat(1:2,:), 2, 1);

rmse_meas = sqrt(mean(err_meas.^2));
rmse_kf   = sqrt(mean(err_kf.^2));

fprintf("Scenario 1 (No dropouts)\n");
fprintf("  TRUE sigma_a: %.2f, KF sigma_a_filt: %.2f\n", sigma_a_true, sigma_a_filt);
fprintf("  TRUE sigma_z: %.2f, KF sigma_z_filt: %.2f\n", sigma_z_true, sigma_z_filt);
fprintf("  RMSE (Measurement): %.3f m\n", rmse_meas);
fprintf("  RMSE (Kalman):      %.3f m\n", rmse_kf);

% Plots
figure; hold on; grid on;
plot(Xtrue(1,:), Xtrue(2,:), 'LineWidth', 2);
plot(xhat(1,:),  xhat(2,:),  'LineWidth', 2);
scatter(Z(1,:), Z(2,:), 15, 'filled');
legend("True", "KF estimate", "GPS meas.", "Location","best");
xlabel("x (m)"); ylabel("y (m)");
title("Scenario 1: 2D Tracking (No Dropouts) - LS Init");

figure; hold on; grid on;
plot(t, err_kf, 'LineWidth', 2);
plot(t, err_meas, 'LineWidth', 1);
legend("KF position error", "GPS position error", "Location","best");
xlabel("Time (s)"); ylabel("Position error magnitude (m)");
title("Scenario 1: Error vs Time (No Dropouts)");

% ---- Save Scenario 1 data for direct comparison in Scenario 2 ----
Z_full = Z;       
xhat1  = xhat;     % scenario 1 estimate
save('scenario1_data.mat', 'Xtrue','Z_full','xhat1','t','dt','N', ...
     'sigma_a_true','sigma_z_true','sigma_a_filt','sigma_z_filt');
