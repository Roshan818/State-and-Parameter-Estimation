function double_pendulum_trajectory
    clc
    clear
    close all
    warning off

    %% Define Constants and Parameters
    g = 9.81; % Acceleration due to gravity (m/s^2)
    m1 = 1;   % Mass of first pendulum (kg)
    m2 = 1;   % Mass of second pendulum (kg)
    L1 = 1;   % Length of first pendulum (m)
    L2 = 1;   % Length of second pendulum (m)

    % Define simulation parameters
    T = 10;       % Total simulation time (s)
    freq = 75;   % Sampling frequency (Hz)
    dt = 1 / freq;
    time = 0:dt:T-dt;
    N = T * freq;

    % Initial conditions
    q1_0 = pi / 2;    % Initial angle of the first pendulum (rad)
    q2_0 = pi / 4;    % Initial angle of the second pendulum (rad)
    omega1_0 = 0;     % Initial angular velocity of the first pendulum (rad/s)
    omega2_0 = 0;     % Initial angular velocity of the second pendulum (rad/s)

    % State vector [q1; q2; omega1; omega2]
    xi_0 = [q1_0; q2_0; omega1_0; omega2_0];

    %% Simulation
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6); 
    [t, xi] = ode45(@(t, xi) double_pendulum_dynamics(t, xi, g, m1, m2, L1, L2), time, xi_0, options);

    % Extract states
    q1 = xi(:, 1);
    q2 = xi(:, 2);

    % Convert angles to Cartesian coordinates
    x1 = L1 * sin(q1);
    y1 = -L1 * cos(q1);
    x2 = x1 + L2 * sin(q2);
    y2 = y1 - L2 * cos(q2);

    %% Visualization
    figure;
    plot3(x2, y2, zeros(size(x2)), 'b', 'LineWidth', 2); % Plot trajectory in 3D
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Trajectory of the Pendulum Bob in 3D');
    grid on;
end

function dxdt = double_pendulum_dynamics(t, xi, g, m1, m2, L1, L2)
    q1 = xi(1);
    q2 = xi(2);
    omega1 = xi(3);
    omega2 = xi(4);

    % Lagrangian
    dq1dt = omega1;
    dq2dt = omega2;

    % Compute accelerations
    domega1dt = (m2 * g * sin(q2) * cos(q1 - q2) - m2 * L1 * omega1^2 * sin(q1 - q2) * cos(q1 - q2) - (m1 + m2) * g * sin(q1)) / (L1 * (m1 + m2 * sin(q1 - q2)^2));
    domega2dt = ((m1 + m2) * (L1 * omega1^2 * sin(q1 - q2) - g * sin(q2) + g * sin(q1) * cos(q1 - q2)) - m2 * L2 * omega2^2 * sin(q1 - q2)) / (L2 * (m1 + m2 * sin(q1 - q2)^2));

    dxdt = [dq1dt; dq2dt; domega1dt; domega2dt];
end
