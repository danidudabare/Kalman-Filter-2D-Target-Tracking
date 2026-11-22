% Scenario 2: 2D Kalman Filter Tracking (With GPS Dropouts)
clear; clc; close all;

% Load Scenario 1 data (for fair comparison)
load('scenario1_data.mat', 'Xtrue','Z_full','xhat1','t','dt','N', ...
     'sigma_a_true','sigma_z_true','sigma_a_filt','sigma_z_filt');

% Dropout settings
dropout_prob = 0.5;  % 50% missing measurements
rng(21);              % separate seed so dropouts are repeatable

% Model
F = [1 0 dt 0;
     0 1 0 dt;
     0 0 1  0;
     0 0 0  1];

H = [1 0 0 0;
     0 1 0 0];

% Q and R (use SAME filter assumptions as Scenario 1)
Q = sigma_a_filt^2 * [dt^4/4   0      dt^3/2  0;
                      0      dt^4/4   0      dt^3/2;
                      dt^3/2  0      dt^2    0;
                      0      dt^3/2   0      dt^2];

R = sigma_z_filt^2 * eye(2);

% Apply GPS dropouts to the SAME measurements used in Scenario 1
Z = Z_full;
dropMask = rand(1, N) < dropout_prob;   % 1 = missing
Z(:, dropMask) = nan;

% Kalman Filter (skip update if missing)
xhat2 = zeros(4, N);
P     = zeros(4,4,N);
I = eye(4);

% ---------- LS init using first M AVAILABLE measurements ----------
M = 10;
avail = find(~isnan(Z(1,:)));     % indices where GPS exists

if numel(avail) < 2
    error("Not enough measurements to initialize.");
end

M_use = min(M, numel(avail));
idxM = avail(1:M_use);

A = [ones(M_use,1), t(idxM)'];    % fit pos = p0 + v*t using absolute time

bx = A \ Z(1,idxM)';              % [x0_est; vx_est] at t=0
by = A \ Z(2,idxM)';              % [y0_est; vy_est] at t=0

xhat2(:,1) = [bx(1); by(1); bx(2); by(2)];
P(:,:,1)   = diag([sigma_z_filt^2, sigma_z_filt^2, 20^2, 20^2]);

% ---------- KF recursion ----------
for k = 2:N
    % Predict
    x_pred = F*xhat2(:,k-1);
    P_pred = F*P(:,:,k-1)*F' + Q;

    % Update only if measurement exists
    if ~isnan(Z(1,k))
        y = Z(:,k) - H*x_pred;
        S = H*P_pred*H' + R;
        K = P_pred*H' / S;

        xhat2(:,k) = x_pred + K*y;
        P(:,:,k)   = (I - K*H) * P_pred;
    else
        xhat2(:,k) = x_pred;
        P(:,:,k)   = P_pred;
    end
end

% Evaluation (RMSE)
err_kf2 = vecnorm(Xtrue(1:2,:) - xhat2(1:2,:), 2, 1);

% measurement error only where Z exists
err_meas2 = nan(1, N);
idx = ~isnan(Z(1,:));
err_meas2(idx) = vecnorm(Xtrue(1:2,idx) - Z(:,idx), 2, 1);

rmse_meas2 = sqrt(nanmean(err_meas2.^2));
rmse_kf2   = sqrt(mean(err_kf2.^2));

fprintf("Scenario 2 (With dropouts)\n");
fprintf("  Dropout probability: %.2f\n", dropout_prob);
fprintf("  TRUE sigma_a: %.2f, KF sigma_a_filt: %.2f\n", sigma_a_true, sigma_a_filt);
fprintf("  TRUE sigma_z: %.2f, KF sigma_z_filt: %.2f\n", sigma_z_true, sigma_z_filt);
fprintf("  RMSE (Measurement, ignoring missing): %.3f m\n", rmse_meas2);
fprintf("  RMSE (Kalman):                        %.3f m\n", rmse_kf2);

% Plots

% Trajectory comparison (Scenario 1 KF vs Scenario 2 KF)
figure; hold on; grid on;
plot(Xtrue(1,:), Xtrue(2,:), 'LineWidth', 2);
plot(xhat1(1,:), xhat1(2,:), 'LineWidth', 2);
plot(xhat2(1,:), xhat2(2,:), 'LineWidth', 2);
scatter(Z(1,idx), Z(2,idx), 15, 'filled'); % only available measurements
legend("True", "KF estimate (Scenario 1: no dropouts)", ...
              "KF estimate (Scenario 2: dropouts)", ...
              "GPS meas. (available)", "Location","best");
xlabel("x (m)"); ylabel("y (m)");
title("Scenario 1 vs Scenario 2: 2D Tracking (Effect of GPS Dropouts)");

% 2) Error vs time
err_kf1 = vecnorm(Xtrue(1:2,:) - xhat1(1:2,:), 2, 1);

figure; hold on; grid on;
plot(t, err_kf1, 'LineWidth', 2);
plot(t, err_kf2, 'LineWidth', 2);
legend("KF error (Scenario 1)", "KF error (Scenario 2)", "Location","best");
xlabel("Time (s)"); ylabel("Position error magnitude (m)");
title("KF Position Error vs Time: Scenario 1 vs Scenario 2");

% 3) Dropout timeline
figure;
stairs(t, dropMask, 'LineWidth', 2); grid on;
ylim([-0.1 1.1]);
xlabel("Time (s)"); ylabel("Dropout (1 = missing)");
title("Scenario 2: GPS Dropout Timeline");
