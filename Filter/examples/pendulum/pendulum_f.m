function [new_state] = pendulum_f(state, omega, w, dt)
%PENDULUM_F Propagate state
%
% Syntax:  [new_state] = pendulum_f(state, omega, w, dt)
%
% Inputs:
%    state - state
%    omega - input
%    w - input noise
%    dt - time between propagation
%
% Outputs:
%    new_state - propagated state
g = 9.81; % gravity constant (m/s^2)
L = 1.3; % wire length (m)
m= 1;

mu = 6;
epsilon = 4;
kappa = 3;

e3 = -[0; 0; 1];
e3_i = state.Rot*e3;
u = state.u;
% % Spherical Pendulum
dot_u = [-u(2)*u(3); u(1)*u(3); 0] + g/L*cross(e3, e3_i) + w(4:6);

% % Pendulum cart 
%dot_u = [-u(1); u(2); 0] + g/L*cross(e3, e3_i) + w(4:6);

% % Duffing Oscillator
% dot_u = -mu * dot(u, u)^(epsilon - 1) * u + 0.5 * kappa * dot(u, u) * u;

% % Double Pendulum


new_state.Rot = state.Rot * so3_exp((state.u + w(1:3))*dt);
new_state.u = state.u + dot_u*dt;
end