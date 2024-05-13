clc
clear
close all
warning off

%% 
syms xi [10,1]; syms th  [3,1]
M = 10; m = 0.9; L = 1; 
k = 50; c = 10;
g = 9.81;
e3 = [0;0;1]; C = [1 0; 0 1; 0 0];

x = xi(1:2); % 2x1 cart disp
v = xi(3:4); % 2x1 cart conjugate velocity
q = xi(5:7); % 3x1 pendulum position wrt cart
om = xi(8:10); % 3x1 conjugate angular frequency

Mat1 = [m*L^2*eye(3)  m*L*so3_wedge(q)*C;
        -m*L*C'*so3_wedge(q)  (M+m)*eye(2)];
Mat2 = [m*L*so3_wedge(C*v)*so3_wedge(om)*q + m*g*L*cross(q,e3) + m*L*cross(so3_wedge(om)*C*v,q); 
        m*L*C'*so3_wedge(om)^2*q + k*x + c*v];

a = -inv(Mat1)*Mat2;

dh_pi = [0; 0; 0];
dh_q = [1; 1; 1];

aSO3f = [so3_wedge(om)];   
bSO3f = so3_wedge(dh_pi);

aTSf = a(1:3);   
bTSf = dh_q;

aEf = [v;a(4:5)]; 
bEf = [0; 0; 1; 1];

[L0_aE, L0_bE, L1_aE, L1_bE, L1L1_bE] = tayL0L1(aEf,bEf,xi(1:4));
[L0_aSO3, L0_bSO3, L1_aSO3, L1_bSO3, L1L1_bSO3] = tayL0L1_S2(aSO3f,bSO3f,xi(5:7));
[L0_aTS, L0_bTS, L1_aTS, L1_bTS, L1L1_bTS] = tayL0L1(aTSf,bTSf,xi(8:10));

aE_v = matlabFunction(aEf,'vars',{([xi])});
L0aE_v = matlabFunction(L0_aE,'vars',{([xi])});
L1aE_v = matlabFunction(L1_aE,'vars',{([xi])});
L0bE_v = matlabFunction(L0_bE,'vars',{([xi])});
L1bE_v = matlabFunction(L1_bE,'vars',{([xi])});
% L1L1bE_v = matlabFunction(L1L1_bE,'vars',{([xi])});

aTS_v = matlabFunction(aTSf,'vars',{([xi])});
% bTS_v = matlabFunction(bTSf,'vars',{([xi])});
L0aTS_v = matlabFunction(L0_aTS,'vars',{([xi])});
L1aTS_v = matlabFunction(L1_aTS,'vars',{([xi])});
L0bTS_v = matlabFunction(L0_bTS,'vars',{([xi])});
L1bTS_v = matlabFunction(L1_bTS,'vars',{([xi])});
% L1L1bTS_v = matlabFunction(L1L1_bTS,'vars',{([xi])});

aSO3_v = matlabFunction(aSO3f,'vars',{([xi])});
L0aSO3_v = matlabFunction(L0_aSO3,'vars',{([xi])});
L1aSO3_v = matlabFunction(L1_aSO3,'vars',{([xi])});
L0bSO3_v = matlabFunction(L0_bSO3,'vars',{([xi])});
L1bSO3_v = matlabFunction(L1_bSO3,'vars',{([xi])});
% L1L1bSO3_v = matlabFunction(L1L1_bSO3,'vars',{([xi])});

%% Parameters
% sequence time (s)
NSIM = 1; % number of ensembles
T = 50; % Time 
dt = 0.01;
sigWE = 0.3; sigZE = 0.2;
sigWSO3 = 0.05; sigZSO3 = 0.05;
sigWTS = 0.01; sigZTS = 0.01;

deltamat = [sqrt(dt)            0;
      dt^1.5/2    dt^1.5/(2*sqrt(3))];
freq = 1/dt;
time = 0:dt:T-dt;
N = T*freq;

for nsim = 1:NSIM
nsim
y = zeros(10,N);

y(:,1) = [0;0; % x
          0;-1; % xdot
          0;0;-1;  % q
         0.5;0.2;0]; % omega

for n = 1:N
  n;
    DWE = sigWE*(deltamat(1,:)*randn(2,4))'; 
    DZE = sigZE*(deltamat(2,:)*randn(2,4))'; 
    DWSO3 = sigWSO3*(deltamat(1,:)*randn(2,3))'; 
    DZSO3 = sigZSO3*(deltamat(2,:)*randn(2,3))'; 
    DWTS = sigWTS*(deltamat(1,:)*randn(2,1))'; 
    DZTS = sigZTS*(deltamat(2,:)*randn(2,1))'; 

% Motion of cart (Euclidean)
    aE = aE_v([y(:,n)]);
    bE = bEf;
    L0aE = L0aE_v([y(:,n)]);
    L1aE = L1aE_v([y(:,n)]);
    L0bE = L0bE_v([y(:,n)]);
    L1bE = L1bE_v([y(:,n)]);
%     L1L1bE = L1L1bE_v([y(:,n)]);

    y(1:4,n+1) = y(1:4,n) + aE.*dt + bE.*DWE...
                      + L0aE.*dt^2/2 + L1bE.*(DWE.^2-dt)/2 + L1aE.*DZE...
                      + L0bE.*(DWE*dt - DZE);% + L1L1bE.*((1/3)*DWE - dt).*DWE;
    y(1:4,n) = y(1:4,n+1);    

% Morion of position vector of pendulum (S2)
    aSO3 = aSO3_v([y(:,n)]);
    bSO3 = bSO3f;
    L0aSO3 = L0aSO3_v([y(:,n)]);
    L1aSO3 = L1aSO3_v([y(:,n)]);
    L0bSO3 = L0bSO3_v([y(:,n)]);
    L1bSO3 = L1bSO3_v([y(:,n)]);
%     L1L1bSO3 = L1L1bSO3_v([y(:,n)]);


    OMEG =  aSO3.*dt + bSO3.*DWSO3...
                  + L0aSO3.*dt^2/2 + L1bSO3.*(DWSO3.^2-dt)/2 + L1aSO3.*DZSO3...
                  + L0bSO3.*(DWSO3*dt - DZSO3);% + L1L1bSO3.*((1/3)*DWSO3 - dt).*DWSO3;

    y(5:7,n+1) = so3_exp_new(OMEG)*y(5:7,n);   

% Morion of angular momentum of pendulum (TSO3)
    aTS = aTS_v([y(:,n)]);
    bTS = bTSf;%bTS_v([y(:,n)]);
    L0aTS = L0aTS_v([y(:,n)]);
    L1aTS = L1aTS_v([y(:,n)]);
    L0bTS = L0bTS_v([y(:,n)]);
    L1bTS = L1bTS_v([y(:,n)]);
%     L1L1bTS = L1L1bTS_v([y(:,n)]);

    y(8:10,n+1) = y(8:10,n) + aTS.*dt + bTS.*DWTS...
                      + L0aTS.*dt^2/2 + L1bTS.*(DWTS.^2-dt)/2 + L1aTS.*DZTS...
                      + L0bTS.*(DWTS*dt - DZTS);% + L1L1bTS.*((1/3)*DWTS - dt).*DWTS;

end
YYc(:,:,nsim) = y(1:2,:);
YYvc(:,:,nsim) = y(3:4,:);
YYp(:,nsim,:) = y(5:7,:);
momp(:,:,nsim) = y(8:10,:);

clear y
end

% Taking the mean of responses
Yc_GTay = nanmean(YYc,3);
Yvc_GTay = nanmean(YYvc,3);
for ii = 1:N+1
Yp_GTay(:,ii) = karcher_mean_sphere(YYp(:,:,ii)); % Karchar's ensemble mean
end
momp_GTay = nanmean(momp,3);


% storing final responses
yC =Yc_GTay; % Position of cart 
Yv = Yvc_GTay; % Velocity of cart 
Q_r =Yp_GTay; % Position of pendulum (Local)
OM = momp_GTay; % Angular velocity of pendulum
yP = C*yC + L*Q_r; % Position of pendulum (Global)

%% Plotting
figure,
for i = 1:50:length(time)
    i;
plot3(yP(1,1:i),yP(2,1:i),yP(3,1:i),'k','linewidth',2); hold on
view(128.5,29); grid on
plot(yC(1,1:i),yC(2,1:i),'r','linewidth',2); 
plot3([yC(1,i);yP(1,i)],[yC(2,i);yP(2,i)],[0;yP(3,i)],'linewidth',2)

[xx yy] = meshgrid(-1.5:0.1:1.5); % Generate x and y data
zz = zeros(size(xx, 1)); % Generate z data
FIG = surf(xx, yy, zz); % Plot the surface
set(FIG, 'FaceAlpha', 0.5,'EdgeAlpha', 0.2)

xlim([-1.5 1.5]); ylim([-1.6 1.6]); zlim([-1.5 1.5]);
drawnow
hold off

LL(i) = sqrt((yC(1,i)-yP(1,i))^2 + (yC(2,i)-yP(2,i))^2 + (0-yP(3,i))^2);
end

%%
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
plot3(Q_r(1,:),Q_r(2,:),Q_r(3,:),'k','linewidth',2)

%
% for i = 2:100:length(time)
%    plot3(q2(1,1:i),q2(2,1:i),q2(3,1:i),'k','linewidth',2);%hold on
%    drawnow
% end
