classdef RK4_Double_pend
methods(Static)
%% Function 1
function dxdt = RK4(ti,xi,h)
%
% ________________________________________________________________________    
k1=RK4_Double_pend.parameters(ti,xi);
k2=RK4_Double_pend.parameters(ti+h/2,xi+h/2*k1);
k3=RK4_Double_pend.parameters(ti+h/2,xi+h/2*k2);
k4=RK4_Double_pend.parameters(ti+h,xi+h*k1);
dxdt=xi+h/6*(k1+2*k2+2*k3+k4);
end


function dxdt = parameters(t,xi)
%
% syms xi [10,1]
m1 = 5; m2 = 1; L1 = 1; L2 = 0.5;
g = 9.81;

    q1 = xi(1:2); q2 = xi(3:5); omeg1 = xi(6); omeg2 = xi(7:9);

%%
    
    S = [0 -1;
         1 0];
    Q = [0 0; 1 0; 0 1];
    e3 = [0;0;1];
     
    mat1 = [(m1+m2)*L1^2   -m2*L1*L2*q1'*S'*Q'*so3_wedge(q2);
            -m2*L1*L2*so3_wedge(q2)'*Q*S*q1  m2*L2^2*eye(3)];
        
    matv2 = -Q*q1*omeg1^1 + omeg1*so3_wedge(q2)*so3_wedge(Q*S*q1)...
            *so3_wedge(omeg2)*q2 + omeg1*so3_wedge(omeg2)*Q*S*q1;
    
    mat2 = [-m2*L1*L2*q1'*S*Q'*so3_wedge(omeg2)^2*q2;
            -cross(m2*L1*L2*q2, matv2)];
        
    mat3 = [(m1+m2)*g*L1*e3'*Q*S*q1;
            -cross(m2*g*L2*q2, e3)];    
        
        
    dxdt=[omeg1*S*q1;
         so3_wedge(omeg2)*q2;
         -mat1\(mat2+mat3)];

end

end
end