clc
clear
close all

%% Initialization
% sequence time (s)
T = 10;
% model frequency (Hz)
freq = 1000;
dt = 1/freq;
time = 0:dt:T-dt;
N = T*freq;

m1 = 5; m2 = 1; L1 = 1; L2 = 0.5;
g = 9.81;  Q = [0 0; 1 0; 0 1];

q = zeros(9,N); 
q(:,1) = [L1;0;0;L2;0;0.5;0;0.2;0.8];

for i = 1:N
    i    
    dxdt = RK4_Double_pend.RK4(time(i),q(:,i),dt);
    q(:,i+1) = dxdt;
%     om_gam_dot(i) = dot(nf(q(:,i)),q(4:6,i)); % check
end

% storing variables
q1 = q(1:2,:); 
q2 = -q(3:5,:); 
omeg1 = q(6,:); 
omeg2 = q(7:9,:);

x1 = L1*Q*q1;
x2 = x1 + L2*q2;
%% Plotting
figure,
plot(q1(1,:),q1(2,:),'k','linewidth',2)

% for i = 2:length(time)
%     plot(q1(1,i),q1(2,i),'ko','linewidth',2,'MarkerSize', 12);%hold on
%     xlim([-1 1]); ylim([-1 1]);
%     drawnow
% end

%%
[X1,Y1,Z1] = sphere(25);
X2 = X1 * L2;
Y2 = Y1 * L2;
Z2 = Z1 * L2;
% 
figure,
FIG = surf(X2,Y2,Z2);
set(FIG, 'FaceAlpha', 0.5,'EdgeAlpha', 0.2)
% shading interp
axis equal
hold on
plot3(q2(1,:),q2(2,:),q2(3,:),'k','linewidth',2)
% 
% for i = 2:100:length(time)
%     plot3(q2(1,1:i),q2(2,1:i),q2(3,1:i),'k','linewidth',2);%hold on
%     drawnow
% end

%% Plotting

figure,
for i = 1:50:length(time)
    i;
plot3(x2(1,1:i),x2(2,1:i),x2(3,1:i),'k','linewidth',2); hold on
view(128.5,29); grid on
plot3(zeros(1,length(q(1,1:i))),q(1,1:i),q(2,1:i),'r','linewidth',4); 
plot3([0;x2(1,i)],[q(1,i);x2(2,i)],[q(2,i);x2(3,i)],'linewidth',2)

% [xx yy] = meshgrid(-1.5:0.1:1.5); % Generate x and y data
% zz = zeros(size(xx, 1)); % Generate z data
% FIG = surf(xx, yy, zz); % Plot the surface
% set(FIG, 'FaceAlpha', 0.5,'EdgeAlpha', 0.2)

xlim([-0.5 0.5]); ylim([-1.6 1.6]); zlim([-1.5 0]);
drawnow
hold off

% LL(i) = sqrt((yC(1,i)-yP(1,i))^2 + (yC(2,i)-yP(2,i))^2 + (0-yP(3,i))^2);
end
%% 
figure,
for i = 1:10:length(time)
plot3(x1(1,:),x1(2,:),x1(3,:),'r','linewidth',0.5); hold on
plot3(x2(1,1:i),x2(2,1:i),x2(3,1:i),'g','linewidth',0.5); hold on
plot3(x1(1,i),x1(2,i),x1(3,i),'ro','linewidth',2); hold on
plot3(x2(1,i),x2(2,i),x2(3,i),'ko','linewidth',2);
plot3([x1(1,i);x2(1,i)],[x1(2,i);x2(2,i)],[x1(3,i);x2(3,i)],'linewidth',2)
plot3([0;x1(1,i)],[0;x1(2,i)],[0;x1(3,i)],'linewidth',2)

LL(i) = sqrt((x1(1,i)-x2(1,i))^2 + (x1(2,i)-x2(2,i))^2 + (x1(3,i)-x2(3,i))^2);
xlim([-0.5 0.5]); ylim([-1.5 1.5]); zlim([-1.5 1.5]);
drawnow
hold off
end

%% GIF plotting
close all
fh = figure(5);
fh.WindowState = 'maximized';
axes('Units', 'normalized', 'Position', [0 0 1 1]);
filename = 'double_pend.gif';
% FIG = surf(X2,Y2,Z2);
% set(FIG, 'FaceAlpha', 0.5,'EdgeAlpha', 0.2)
% shading interp
% axis equal
% hold on
for i = 1:100:length(time)
    
plot3(x1(1,:),x1(2,:),x1(3,:),'r','linewidth',1.5); hold on
xlim([-0.5 0.5]); ylim([-1.5 1.5]); zlim([-1.5 1.5]);
set(gca,'visible','off')
plot3(x2(1,1:i),x2(2,1:i),x2(3,1:i),'g','linewidth',0.5);
plot3(x1(1,i),x1(2,i),x1(3,i),'ro','linewidth',8);
plot3(x2(1,i),x2(2,i),x2(3,i),'ko','linewidth',8);
plot3([x1(1,i);x2(1,i)],[x1(2,i);x2(2,i)],[x1(3,i);x2(3,i)],'linewidth',3)
plot3([0;x1(1,i)],[0;x1(2,i)],[0;x1(3,i)],'linewidth',3)

xlim([-0.5 0.5]); ylim([-1.5 1.5]); zlim([-1.5 1.5]);
drawnow
hold off

xlabel('')
set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[]);
set(gcf, 'Color', 'white')
set(gca,'visible','off')
    hold all
%     pause( 0.00001 );
    frame = getframe(5);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1 
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
hold off

% drawnow
% if i == length(ys)
% else
% delete(hh2)
% end
end


