clc
clear
close all
warning off


%% Duffing Oscillator Unscented Kalman Filter
syms xi [6,1]
dof = 3;
mu = 0.9; e3 = [1;0;0]; eps = 1.3; alpha = 5;

q = xi(1:3); % 3x1 pendulum position
pii = xi(4:6); % 3x1 conjugate angular frequency
d_rim = acos(q.'*e3);
logq_e3 = (eye(3)-q*q.')*e3*(d_rim/sin(d_rim));
a1= 1; a2 =1;

H = (1/2)*pii.'*pii + (a1*(1/2)*d_rim^2 + a2*(1/4)*alpha*d_rim^4);

dh_q = [0.0; 0.0; 0.0];
dh_pi = [0.; 0.; 0.];

aSO3f = so3_wedge(pii);   
bSO3f = so3_wedge(dh_q);

aTSf1 = cross(q,(a1+a2*alpha*d_rim^2)*logq_e3)-mu*(cross(pii,q).'*cross(pii,q))^(eps-1)*cross(pii,q);   
bTSf = -cross(q,dh_q)-cross(pii,dh_pi);
aTSf1 = aTSf1 + 0.5*cross(cross(q,bTSf),bTSf);
aTSf = aTSf1-((q.'*aTSf1)/(q.'*q))*q;

[L0_aSO3, L0_bSO3, L1_aSO3, L1_bSO3, L1L1_bSO3] = tayL0L1_S2(aSO3f,bSO3f,xi(1:3));
[L0_aTS, L0_bTS, L1_aTS, L1_bTS, L1L1_bTS] = tayL0L1(aTSf,bTSf,xi(4:6));

var = [xi];
aTS_v = matlabFunction(aTSf,'vars',{(var)});
bTS_v = matlabFunction(bTSf,'vars',{(var)});
L0aTS_v = matlabFunction(L0_aTS,'vars',{(var)});
L1aTS_v = matlabFunction(L1_aTS,'vars',{(var)});
L0bTS_v = matlabFunction(L0_bTS,'vars',{(var)});
L1bTS_v = matlabFunction(L1_bTS,'vars',{(var)});

aSO3_v = matlabFunction(aSO3f,'vars',{(var)});
L0aSO3_v = matlabFunction(L0_aSO3,'vars',{(var)});
L1aSO3_v = matlabFunction(L1_aSO3,'vars',{(var)});
L0bSO3_v = matlabFunction(L0_bSO3,'vars',{(var)});
L1bSO3_v = matlabFunction(L1_bSO3,'vars',{(var)});

%% Parameters
T = 40;
freq = 100;
dt = 1/freq;
time = 0:dt:T-dt;
N = T*freq;
NSIM = 10;
deltamat = [sqrt(dt) 0; dt^1.5/2 dt^1.5/(2*sqrt(3))];
sigWSO3 = 0.5; sigZSO3 = 0.5;
sigWTS = 0.5; sigZTS = 0.5;

% Initialize state and covariance matrices
x_hat = zeros(6, N);
P_hat = eye(6);

for nsim = 1:NSIM
    disp(['Simulation ', num2str(nsim)])
    
    % Initialize state vector and covariance matrix for each simulation
    xi_hat = [0; 0; 1; 0.4; 0; 0]; % Initial state
    P_xi = eye(6);
    
    % Unscented Kalman Filter
    for n = 2:N
        disp(['Time step ', num2str(n)])
        
        % Generate sigma points
        sigma_points = unscented_transform(xi_hat, P_xi);
        
        % Predict
        [xi_hat_minus, P_xi_minus] = ukf_predict(xi_hat, P_xi, dt, sigma_points);
        
        % Update
        y_meas = [q; pii]; % Measurement vector (same as state vector)
        [xi_hat, P_xi] = ukf_update(xi_hat_minus, P_xi_minus, y_meas, dt, sigma_points);
        
        % Store results
        x_hat(:, n) = xi_hat;
    end
end

%% Plotting
L = 1;
[X1,Y1,Z1] = sphere(25);
X2 = X1 * L;
Y2 = Y1 * L;
Z2 = Z1 * L;

figure,
FIG = surf(X2, Y2, Z2);
set(FIG, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.2)
axis equal
hold on
plot3(x_hat(1, 2:end), x_hat(2, 2:end), x_hat(3, 2:end), 'k', 'linewidth', 2)
xlabel('X-axis')
ylabel('Y-axis')
zlabel('Z-axis')
title('Unscented Kalman Filter for Duffing Oscillator')


function sigma_points = unscented_transform(x, P)
    n = length(x);
    alpha = 1e-3;
    beta = 2;
    kappa = 3 - n;
    
    lambda = alpha^2 * (n + kappa) - n;
    
    % Calculate sigma points
    sigma_points(:, 1) = x;
    sqrt_P = sqrtm((n + lambda) * P);
    for i = 1:n
        sigma_points(:, i + 1) = x + sqrt_P(:, i);
        sigma_points(:, i + n + 1) = x - sqrt_P(:, i);
    end
    
    % Adjust weights
    w_m = [lambda / (n + lambda), 0.5 / (n + lambda) * ones(1, 2 * n)];
    w_c = w_m;
    w_c(1) = w_c(1) + (1 - alpha^2 + beta);
    
    % Apply weights to sigma points
    sigma_points = sigma_points * diag(w_c);
end


function [predicted_state, predicted_P] = ukf_predict(state, P, f, phi, Q, weights)
%UKF_PREDICT Prediction step of the UKF
%
% Syntax: [predicted_state, predicted_P] = ukf_predict(state, P, f, phi, Q, weights)
%
% Inputs:
%    state - state
%    P - covariance matrix
%    f - state transition function, function
%    phi - retraction, function
%    Q - process noise covariance matrix
%    weights - weight parameters
%
% Outputs:
%    predicted_state - predicted state
%    predicted_P - predicted covariance matrix

TOL = 1e-9; % tolerance for assuring positivness of P

% set variable size
d = length(P);

P = P + TOL * eye(d);

% set sigma points
w_u = weights.u;
xis = w_u.sqrt_d_lambda * chol(P)';

% propagate sigma points through state transition function
xs = zeros(d, 2*d);
x_hat = f(state);
for j = 1:d
    chi_j_plus = phi(state, xis(:, j));
    chi_j_minus = phi(state, -xis(:, j));
    xs(:, j) = f(chi_j_plus);
    xs(:, d + j) = f(chi_j_minus);
end

% state mean
x_bar = w_u.wm0 * x_hat + w_u.wj * sum(xs, 2);

% prune mean before computing covariance
xs = xs - x_bar;
x_hat = x_hat - x_bar;

% compute covariance matrix
P_xx = w_u.wc0 * (x_hat * x_hat') + w_u.wj * (xs * xs') + Q;

% update state
predicted_state = phi(state, w_u.wj * xs);

% update covariance
predicted_P = P_xx;
predicted_P = (predicted_P + predicted_P') / 2;
end
