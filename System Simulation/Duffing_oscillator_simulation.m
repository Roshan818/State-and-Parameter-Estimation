clc
clear
close all
warning off

%% 
clc
clear
close all
warning off

%% 
syms xi [6,1]
dof = 3;
mu = 0.9; e3 = [1;0;0]; eps = 1.3; alpha = 5;

q = xi(1:3); % 3x1 pendulum position
pii = xi(4:6); % 3x1 conjugate angular frequency
d_rim = acos(q.'*e3);
logq_e3 = (eye(3)-q*q.')*e3*(d_rim/sin(d_rim));
a1= 1; a2 =1;

H = (1/2)*pii.'*pii + (a1*(1/2)*d_rim^2 + a2*(1/4)*alpha*d_rim^4);
%
dh_q = [0.0; 0.0; 0.0];
dh_pi = [0.; 0.; 0.];

aSO3f = so3_wedge(pii);   
bSO3f = so3_wedge(dh_q);
% aSO3f = aSO3f;

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
% L1L1bTS_v = matlabFunction(L1L1_bTS,'vars',{(var)});

aSO3_v = matlabFunction(aSO3f,'vars',{(var)});
L0aSO3_v = matlabFunction(L0_aSO3,'vars',{(var)});
L1aSO3_v = matlabFunction(L1_aSO3,'vars',{(var)});
L0bSO3_v = matlabFunction(L0_bSO3,'vars',{(var)});
L1bSO3_v = matlabFunction(L1_bSO3,'vars',{(var)});
% L1L1bSO3_v = matlabFunction(L1L1_bSO3,'vars',{(var)});

%% Parameters
% sequence time (s)
T = 40;
freq = 100;
dt = 1/freq;
time = 0:dt:T-dt;
N = T*freq;
NSIM = 10;
deltamat = [sqrt(dt)            0;
          dt^1.5/2    dt^1.5/(2*sqrt(3))];
sigWSO3 = 0.5; sigZSO3 = 0.5;
sigWTS = 0.5; sigZTS = 0.5;
for nsim = 1:NSIM
    nsim
    y = zeros(6,N);

        y(:,2) = [0;0;1;  % q
                 0.4;0;0]; % omega
     %%%%%    Simulation over samples
    for n = 2:N
      n;
        DWSO3 = sigWSO3*(deltamat(1,:)*randn(2,3))'; 
        DZSO3 = sigZSO3*(deltamat(2,:)*randn(2,3))'; 
        DWTS = sigWTS*(deltamat(1,:)*randn(2,1))'; 
        DZTS = sigZTS*(deltamat(2,:)*randn(2,1))'; 

    % Morion of position vector of pendulum (SO3)
    varsS2 = [y(1:3,n);y(4:6,n)];
        aSO3 = aSO3_v(varsS2);
        bSO3 = bSO3f;
        L0aSO3 = L0aSO3_v(varsS2);
        L1aSO3 = L1aSO3_v(varsS2);
        L0bSO3 = L0bSO3_v(varsS2);
        L1bSO3 = L1bSO3_v(varsS2);
    
OMEG =  aSO3.*dt + bSO3.*DWSO3...
       + L0aSO3.*dt^2/2 + L1bSO3.*(DWSO3.^2-dt)/2 + L1aSO3.*DZSO3...
       + L0bSO3.*(DWSO3*dt - DZSO3);
y(1:3,n+1) = so3_exp_new(OMEG)*y(1:3,n);
        
    
    % Morion of angular momentum of pendulum (TSO3)
        varsTS = [y(1:3,n+1);y(4:6,n)];
        aTS = aTS_v(varsTS);
        bTS = bTS_v(varsTS);
        L0aTS = L0aTS_v(varsTS);
        L1aTS = L1aTS_v(varsTS);
        L0bTS = L0bTS_v(varsTS);
        L1bTS = L1bTS_v(varsTS);

 omg_upd = y(4:6,n) + (aTS.*dt + bTS.*DWTS...
              + L0aTS.*dt^2/2 + L1bTS.*(DWTS.^2-dt)/2 + L1aTS.*DZTS...
              + L0bTS.*(DWTS*dt - DZTS));
y(4:6,n+1) = ParTransport(y(1:3,n),y(1:3,n+1))*omg_upd;
    end
%

Q_rn(:,nsim,:) = y(1:3,:);
PIn(:,:,nsim) = y(4:6,:);
end
% clearvars -except Q_rn PIn Ham_exp_n Ham_n m L N time dt 

%% storing variables
for ii = 1:N+1
Q_r(:,ii) = karcher_mean_sphere(Q_rn(:,:,ii)); % Karchar's ensemble mean
end
PI = nanmean(PIn,3);
OM = PI;

%%
L = 1;
[X1,Y1,Z1] = sphere(25);
X2 = X1 * L;
Y2 = Y1 * L;
Z2 = Z1 * L;
%
figure,
FIG = surf(X2,Y2,Z2);
set(FIG, 'FaceAlpha', 0.5,'EdgeAlpha', 0.2)
% shading interp
axis equal
hold on
plot3(Q_r(1,2:end),Q_r(2,2:end),Q_r(3,2:end),'k','linewidth',2)
