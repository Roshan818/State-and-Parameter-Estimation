%% Initialization
% Start by cleaning the workspace.
clc;
clear;
close all;

%% Model and Simulation
% sequence time (s)
T = 10;
% model frequency (Hz)
freq = 100;
% model noise standard deviation (noise is isotropic)
model_noise_std = [1/180*pi;  % orientation (rad) 
                  1/180*pi];  % orientation velocity (rad/s)
%% simulate true states and noisy inputs
% total number of timestamps
N = T*freq;
% integration step (s)
dt = 1/freq;
% set noise to zero to compute the true trajectory
w = zeros(6, 1);
rpy = [57.3/180*pi; 40/180*pi; 0];

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

%% Plotting
[X,Y,Z] = sphere(25);
r = L;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;
% 
figure(1),
FIG = surf(X2,Y2,Z2);
set(FIG, 'FaceAlpha', 0.5,'EdgeAlpha', 0.2)
% shading interp
axis equal
hold on
plot3(ys(1,:),ys(2,:),ys(3,:),'k','linewidth',2); % Spherical 
plot(ys(1,:),ys(2,:),'k','linewidth',2) 


for i = 1:length(ys)
 plot(ys(1,1:i),ys(2,1:i),'k','linewidth',2)
 drawnow;
end

