function [L0a, L0b, L1a, L1b, L1L1b] = tayL0L1_S2(a,b,y)
syms t
La = size(b,1); Lb = size(b,2);
bsq = b^2;

%% L0a
    misc = sym(zeros(size(a)));
for i = 1:La
    for j = 1:La
    misc = misc + a(i,j)*diff(a,y(j))*y(j) +0.5*bsq(i,j)*y(j)*diff(a,y(j));
    end
end
for k = 1:Lb
    for l = 1:Lb
        for i = 1:La
            for j = 1:La
                misc = misc+0.5*b(i,j)*y(j)*b(k,l)*y(l)...
                      *diff(diff(a,y(j)),y(l));
            end
        end
    end
end
    L0a = misc;
fprintf("L0a is created.\n")
%% L0b
    misc = sym(zeros(size(b)));
for i = 1:La
    for j = 1:La
    misc = misc + a(i,j)*diff(b,y(j))*y(j) +0.5*bsq(i,j)*y(j)*diff(b,y(j));
    end
end
for k = 1:Lb
    for l = 1:Lb
        for i = 1:La
            for j = 1:La
                misc = misc+0.5*b(i,j)*y(j)*b(k,l)*y(l)...
                      *diff(diff(b,y(j)),y(l));
            end
        end
    end
end
    L0b = misc;
fprintf("L0b is created.\n")

%% Lja
    misc = sym(zeros(size(a)));  
for i = 1:La
    for j = 1:La
        misc = misc + b(i,j)*diff(a,y(j))*y(j);
    end
end
L1a = misc;
fprintf("L1a is created.\n")

%% Ljb
    misc = sym(zeros(size(b)));
for i = 1:La
    for j = 1:La
        misc = misc + b(i,j)*diff(b,y(j))*y(j);
    end
end
L1b = misc;
fprintf("L1b is created.\n")

%% L1L1b
    misc = sym(zeros(size(b)));
for i = 1:La
    for j = 1:La
        misc = misc + b(i,j)*diff(L1b,y(j))*y(j);
    end
end
L1L1b = misc;
fprintf("L1L1b is created.\n")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% L0a
%     misc = sym(zeros(size(a)));
% for i = 1:La
%     for j = 1:La
%     misc = misc + a(i,j)*diff(a,y(j))*y(j);
%     end
% end
% for k = 1:Lb
%     for l = 1:Lb
%         for i = 1:La
%             for j = 1:La
%                 misc = misc+0.5*b(i,j)*b(k,l)...
%                       *diff(diff(a,y(j))*y(j),y(l))*y(l);
%             end
%         end
%     end
% end
%     L0a = misc;
% fprintf("L0a is created.\n")
% 
% %% L0b
%     misc = sym(zeros(size(b)));
% for i = 1:La
%     for j = 1:La
%     misc = misc + a(i,j)*diff(b,y(j))*y(j);
%     end
% end
% for k = 1:Lb
%     for l = 1:Lb
%         for i = 1:La
%             for j = 1:La
%                 misc = misc+0.5*b(i,j)*b(k,l)...
%                       *diff(diff(b,y(j))*y(j),y(l))*y(l);
%             end
%         end
%     end
% end
%     L0b = misc;
% fprintf("L0b is created.\n")
% 
% %% Lja
%     misc = sym(zeros(size(a)));  
% for i = 1:La
%     for j = 1:La
%         misc = misc + b(i,j)*diff(a,y(j))*y(j);
%     end
% end
% L1a = misc;
% fprintf("L1a is created.\n")
% 
% %% Ljb
%     misc = sym(zeros(size(b)));
% for i = 1:La
%     for j = 1:La
%         misc = misc + b(i,j)*diff(b,y(j))*y(j);
%     end
% end
% L1b = misc;
% fprintf("L1b is created.\n")
% 
% %% L1L1b
%     misc = sym(zeros(size(b)));
% for i = 1:La
%     for j = 1:La
%         misc = misc + b(i,j)*diff(L1b,y(j))*y(j);
%     end
% end
% L1L1b = misc;
% fprintf("L1L1b is created.\n")
end

