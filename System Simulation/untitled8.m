clc;
clear;
close all;

%% Define Parameters and Initial Conditions
T = 10;             % Total time
freq = 1000;        % Model frequency
dt = 1/freq;        % Time step
N = T*freq;         % Number of time steps

% System parameters
m1 = 5; m2 = 1;     % Masses
L1 = 1; L2 = 0.5;   % Lengths
g = 9.81;           % Gravity

% Initial conditions [q1_x, q1_y, q2_x, q2_y, theta1, theta2, omega1, omega2, omega3]
q0 = [L1; 0; L2; 0; 0; 0; 0.5; 0.2; 0.8]; 

% Measurement noise covariance matrix
R = 0.02^2 * eye(2);

%% Model Dynamics
% Define the dynamics equations for the double pendulum system
% This function describes how the state variables change over time
f_double_pendulum = @(q) double_pendulum_dynamics(q, m1, m2, L1, L2, g);

%% Measurement Function (h)
% Measurement function maps the state variables to the measurements
% In this case, we directly observe the positions of the pendulum segments
h_double_pendulum = @(q) q(1:4);

%% UKF Initialization
% Initial uncertainty matrix
P0 = diag([1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3]);

% Define UKF parameters
alpha = 1e-3;
beta = 2;
kappa = 0;
weights = ukf_set_weight(9, 2, alpha, beta, kappa);

% Initialize state estimate and covariance
ukf_state = q0;
ukf_P = P0;

% Storage for estimated states and covariances
ukf_states = zeros(length(q0), N);
ukf_Ps = zeros(length(q0), length(q0), N);

%% Filtering
for n = 1:N
    % Prediction Step
    [ukf_state, ukf_P] = ukf_propagation(f_double_pendulum, ukf_state, ukf_P, dt, weights);
    
    % Update Step (No measurements for this example)
    % Skipping measurement update as we don't have measurements in this example
    
    % Store estimates
    ukf_states(:, n) = ukf_state;
    ukf_Ps(:, :, n) = ukf_P;
end

%% Plotting Results
time = (0:dt:T-dt)';
figure;
plot(time, ukf_states(1,:), 'b', 'LineWidth', 2);
hold on;
plot(time, ukf_states(3,:), 'r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Position');
legend('Segment 1', 'Segment 2');
title('Estimated Positions of Pendulum Segments');
grid on;

%% Dynamics Function
% Dynamics function for the double pendulum system
function dqdt = double_pendulum_dynamics(q, m1, m2, L1, L2, g)
    % Unpack state variables
    q1_x = q(1); q1_y = q(2); q2_x = q(3); q2_y = q(4);
    theta1 = q(5); theta2 = q(6);
    omega1 = q(7); omega2 = q(8); omega3 = q(9);
    
    % Dynamics equations
    dqdt = zeros(9, 1);
    dqdt(1) = omega1;
    dqdt(2) = omega2;
    dqdt(3) = omega3;
    dqdt(4) = - (m1+m2) * g / m1 * q1_y - m2 * g / m1 * q2_y;
    dqdt(5) = omega1;
    dqdt(6) = omega2;
    dqdt(7) = 0; % Placeholder, needs to be updated with proper dynamics
    dqdt(8) = 0; % Placeholder, needs to be updated with proper dynamics
    dqdt(9) = 0; % Placeholder, needs to be updated with proper dynamics
end

%% UKF Propagation Function
function [state_new, P_new] = ukf_propagation(f, state, P, dt, weights)
    % Number of state variables
    n = length(state);
    
    % Generate sigma points
    sigma_points = ukf_generate_sigma_points(state, P, weights);
    
    % Propagate sigma points through dynamics function
    sigma_points_pred = zeros(size(sigma_points));
    for i = 1:size(sigma_points, 2)
        sigma_points_pred(:,i) = f(sigma_points(:,i));
    end
    
    % Predicted state and covariance
    state_pred = sum(bsxfun(@times, sigma_points_pred, weights.wm), 2);
    P_pred = zeros(n); % Correct initialization
    for i = 1:size(sigma_points_pred, 2)
        P_pred = P_pred + weights.cov(i) * (sigma_points_pred(:,i) - state_pred) * (sigma_points_pred(:,i) - state_pred)';
    end
    
    % Add process noise covariance (Q) here if applicable
    
    % Output
    state_new = state_pred;
    P_new = P_pred;
end


%% Helper Functions
function sigma_points = ukf_generate_sigma_points(state, P, weights)
    n = length(state);
    sigma_points = zeros(n, 2*n+1);
    sigma_points(:, 1) = state;
    sigma_sqrt = sqrtm(P);
    for i = 1:n
        sigma_points(:, i+1) = state + sigma_sqrt(:, i) * weights.sqrt_lambda;
        sigma_points(:, i+n+1) = state - sigma_sqrt(:, i) * weights.sqrt_lambda;
    end
end

function weights = ukf_set_weight(n, m, alpha, beta, kappa)
    lambda = alpha^2 * (n + kappa) - n;
    c = n + lambda;
    wm = zeros(1, 2*n+1);
    wc = zeros(1, 2*n+1);
    wm(1) = lambda / c;
    wc(1) = lambda / c + (1 - alpha^2 + beta);
    for i = 2:2*n+1
        wm(i) = 1 / (2*c);
        wc(i) = 1 / (2*c);
    end
    weights.wm = wm;
    weights.wc = wc;
    weights.sqrt_lambda = sqrt(c);
    weights.cov = [wm(1) * ones(1, n), (wm(1) + (1 - alpha^2 + beta) * c) * ones(1, n)];
end
