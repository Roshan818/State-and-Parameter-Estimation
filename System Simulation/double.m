%% Initialization
clc;
clear;
close all;

%% Model and Simulation Parameters
T = 10; % Total simulation time (s)
freq = 100; % Model frequency (Hz)
model_noise_std = [1/180*pi;  % Orientation (rad) 
                   1/180*pi;  % Orientation velocity (rad/s)
                   1/180*pi;  % Second pendulum orientation (rad)
                   1/180*pi]; % Second pendulum orientation velocity (rad/s)

%% Simulate True States and Noisy Inputs
N = T * freq; % Total number of timestamps
dt = 1 / freq; % Integration step (s)
w = zeros(10, 1); % Noise vector initialization

% Initial conditions
rpy = [57.3/180*pi; 40/180*pi; 0; 20/180*pi; 0]; % Initial roll, pitch, yaw angles (rad)
states(1).Rot1 = so3_from_rpy(rpy(1:3)); % Initial rotation matrix for first pendulum
states(1).Rot2 = so3_from_rpy(rpy(4:5)); % Initial rotation matrix for second pendulum
states(1).u1 = [-10/180*pi; 30/180*pi; 0]; % Initial angular velocities for first pendulum (rad/s)
states(1).u2 = [-10/180*pi; 30/180*pi; 0]; % Initial angular velocities for second pendulum (rad/s)

% Simulation loop
for n = 2:N
    % Generate noise
    w(1:5) = model_noise_std .* randn(5, 1);
    
    % System parameters
    g = 9.81; % Gravity constant (m/s^2)
    L1 = 1.3; % Length of first pendulum (m)
    L2 = 1.3; % Length of second pendulum (m)
    
    % Compute next state for the first pendulum
    e3 = -[0; 0; 1];
    e3_i = states(n-1).Rot1 * e3;
    u1 = states(n-1).u1;
    dot_u1 = [-u1(2)*u1(3); u1(1)*u1(3); 0] + g/L1 * cross(e3, e3_i) + w(1:3);
    states(n).Rot1 = states(n-1).Rot1 * so3_exp((states(n-1).u1 + w(1:3)) * dt);
    states(n).u1 = states(n-1).u1 + dot_u1 * dt;
    
    % Compute next state for the second pendulum
    e3_i2 = states(n-1).Rot2 * e3;
    u2 = states(n-1).u2;
    dot_u2 = [-u2(2)*u2(3); u2(1)*u2(3); 0] + g/L2 * cross(e3, e3_i2) + w(4:5); % Adjusted indexing
    states(n).Rot2 = states(n-1).Rot2 * so3_exp((states(n-1).u2 + w(4:5)) * dt); % Adjusted indexing
    states(n).u2 = states(n-1).u2 + dot_u2 * dt;
    
    % Compute pendulum positions
    ys1(:, n) = L1 * (states(n).Rot1) * e3;
    ys2(:, n) = L1 * (states(n).Rot1) * e3 + L2 * (states(n).Rot2) * e3;
end

%% Plotting
[X, Y, Z] = sphere(25);
r = L1 + L2;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;

% Plot sphere
figure(1);
FIG = surf(X2, Y2, Z2);
set(FIG, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.2);
axis equal;
hold on;

% Plot pendulum trajectory
plot3(ys1(1,:), ys1(2,:), ys1(3,:), 'b', 'linewidth', 2); % Trajectory of first pendulum
plot3(ys2(1,:), ys2(2,:), ys2(3,:), 'r', 'linewidth', 2); % Trajectory of second pendulum

% Animate pendulum trajectories
for i = 1:length(ys1)
    plot3(ys1(1,1:i), ys1(2,1:i), ys1(3,1:i), 'b', 'linewidth', 2);
    plot3(ys2(1,1:i), ys2(2,1:i), ys2(3,1:i), 'r', 'linewidth', 2);
    drawnow;
end
