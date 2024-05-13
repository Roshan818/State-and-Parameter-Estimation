function [] = pendulum_results_plot(ukf_states, ukf_Ps, states, dt)
%PENDULUM_RESULTS_PLOT
%
% Syntax: [] = pendulum_results_plot(ukf_states, ukf_Ps, states, dt)

set(groot,'defaulttextinterprete','latex');  
set(groot, 'defaultAxesTickLabelInterprete','latex');  
set(groot, 'defaultLegendInterprete','latex');


[Rots, ~] = pendulum_get_states(states);
[ukf_Rots, ~] = pendulum_get_states(ukf_states);

% get position and confidence interval
N = length(states);
t = linspace(0, dt*N, N);
e3 = [0; 0; -1];
e3wedge = so3_wedge(e3);
ps = zeros(3, N);
hat_ps = zeros(3, N);
ukf3sigma = zeros(3, N);
A = eye(6);
L = 1.3;
for n = 1:N
    ps(:, n) = L*Rots{n}*e3;
    hat_ps(:, n) = L*ukf_Rots{n}*e3;
    A(1:3, 1:3) = ukf_Rots{n}*e3wedge;
    P = A*squeeze(ukf_Ps(n, :, :))*A';
    ukf3sigma(:, n) = diag(P(1:3, 1:3));
end
ukf_err = vecnorm(ps-hat_ps);

fig1 = figure;
title('Position')
xlabel('$t$ (s)')
ylabel('position (m)')
hold on;
grid on;
plot(t, ps(1, :),'LineWidth',2);
plot(t, ps(2, :),'LineWidth',2);
plot(t, ps(3, :),'LineWidth',2);
plot(t, hat_ps(1, :),'LineWidth',2);
plot(t, hat_ps(2, :),'LineWidth',2);
plot(t, hat_ps(3, :),'LineWidth',2);
legend('$x$', '$y$', '$z$', '$x$ UKF', '$y$ UKF', '$z$ UKF');
axis equal;

fig2 = figure;
title('Position in $xy$-plan')
xlabel('$x$ (m)')
ylabel('$y$ (m)')
hold on;
grid on;
plot(ps(1, :), ps(2, :), 'r','LineWidth',2);
plot(hat_ps(1, :), hat_ps(2, :), 'b','LineWidth',2);
legend('true position', 'UKF');

fig3 = figure;
title('Position in $yz$-plan')
xlabel('$y$ (m)')
ylabel('$z$ (m)')
hold on;
grid on;
plot(ps(2, :), ps(3, :), 'r','LineWidth',2);
plot(hat_ps(2, :), hat_ps(3, :), 'b','LineWidth',2);
legend('true position', 'UKF');

fig4 = figure;
title('Position in $zx$-plan')
xlabel('$z$ (m)')
ylabel('$x$ (m)')
hold on;
grid on;
plot(ps(3, :), ps(1, :), 'r','LineWidth',2);
plot(hat_ps(3, :), hat_ps(1, :), 'b','LineWidth',2);
legend('true position', 'UKF');


fig5 = figure;
title('Position error (m)')
xlabel('$t$ (s)')
ylabel('position error (m)')
hold on;
grid on;
plot(t, ukf_err, 'b','LineWidth',2);
plot(t, ukf3sigma,'r--','LineWidth',2);
legend('UKF', '$3\sigma$ UKF');

end
