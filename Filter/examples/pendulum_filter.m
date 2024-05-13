%% Pendulum Example

%% Initialization
% Start by cleaning the workspace.
clc
clear

%% Model and Simulation
% sequence time (s)
T = 10;
% model frequency (Hz)
freq = 100;
% model noise standard deviation (noise is isotropic)
model_noise_std = [1/180*pi;  % orientation (rad) 
                  1/180*pi];  % orientation velocity (rad/s)
% set noise to zero to compute the true trajectory
w = zeros(6, 1);
rpy = [57.3/180*pi; 40/180*pi; 0];

% total number of timestamps
N = T*freq;
% integration step (s)
dt = 1/freq;

% init variables at zero and do for loop
omegas(N) = struct;
states(N) = struct;
states(1).Rot = so3_from_rpy(rpy);
states(1).u = [-10/180*pi; 30/180*pi; 0];

for n = 2:N
    w(1:3) = model_noise_std(1)*randn(3, 1);
    w(4:6) = model_noise_std(2)*randn(3, 1);
        
    g = 9.81; % gravity constant (m/s^2)
    L = 1.3; % wire length (m)

    e3 = -[0; 0; 1];
    e3_i = states(n-1).Rot*e3;
    u = states(n-1).u;
    dot_u = [-u(2)*u(3); u(1)*u(3); 0] + g/L*cross(e3, e3_i) + w(4:6);
    states(n).Rot = states(n-1).Rot * so3_exp((states(n-1).u + w(1:3))*dt);
    states(n).u = states(n-1).u + dot_u*dt;

    ys(:, n) = L*(states(n).Rot)*e3;
end

% observation frequency (Hz)
obs_freq = 20;
% observation noise standard deviation (m)
obs_std = 0.02;
% vector to know where measurements happen
one_hot_ys = zeros(N, 1);
one_hot_ys(1:freq/obs_freq:end) = 1; % freq/obs_freq must be integer

%% Filtering
% propagation noise covariance matrix
Q = blkdiag(model_noise_std(1)^2*eye(3), model_noise_std(2)^2*eye(3));
% measurement noise covariance matrix
R = obs_std^2*eye(3);
% initial uncertainty matrix
P0 = blkdiag((45/180*pi)^2*eye(3), (10/180*pi)^2*eye(3));
% sigma point parameters
alpha = [1e-3; 1e-3; 1e-3];
% define the UKF propagation and measurement functions
f = @pendulum_f;
h = @pendulum_h;
phi = @pendulum_phi;
phi_inv = @pendulum_phi_inv;
% get UKF weight parameters
weights = ukf_set_weight(6, 6, alpha);
% compute Cholewski decomposition of Q only once
cholQ = chol(Q);

ukf_state = states(1);
ukf_state.Rot = eye(3);
ukf_state.u = zeros(3, 1);
ukf_P = P0;

% set variables for recording estimates along the full trajectory
ukf_states = ukf_state;
ukf_Ps = zeros(N, length(ukf_P), length(ukf_P));
ukf_Ps(1, :, :) = ukf_P;

%% Filtering
k = 2;
for n = 2:N
    % propagation
    [ukf_state, ukf_P] = ukf_propagation(ukf_state, ukf_P, omegas(n-1), ...
        f, dt, phi, phi_inv, cholQ, weights);
    % update only if a measurement is received
    if one_hot_ys(n) == 1
        [ukf_state, ukf_P] = ukf_update(ukf_state, ukf_P, ys(:, k), ...
            h, phi, R, weights);
        k = k + 1;
    end
    % save estimates
    ukf_states(n) = ukf_state;
    ukf_Ps(n, :, :) = ukf_P;
end

%% Plotting
pendulum_results_plot(ukf_states, ukf_Ps, states, dt);