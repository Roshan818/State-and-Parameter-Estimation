function [L0a, L0b, L1a, L1b, L1L1b] = tayL0L1(a,b,y)

syms t
La = size(b,1); Lb = size(b,2);

%% L0a
    misc = sym(zeros(size(a)));
for j = 1:La
    misc = misc + a(j)*diff(a,y(j));
end
for k = 1:Lb
    for i = 1:La
        for j = 1:La
            misc = misc+0.5*b(i,k)*b(j,k)*diff((diff(a,y(i))),y(j));
        end
    end
end
    L0a = misc;
fprintf("L0a is created.\n")

%% L0b
    misc = sym(zeros(size(b)));
    for i = 1:La
        misc = misc+a(i)*diff(b,y(i));
    end
    for k = 1:Lb
        for i = 1:La
            for j = 1:La
                misc = misc+0.5*b(i,k)*b(j,k)*diff((diff(b,y(i))),y(j));
            end
        end
    end
    L0b = misc;
fprintf("L0b is created.\n")

%% Lja
    misc = sym(zeros(size(b)));
for i = 1:Lb    
    for k = 1:La
        for j = 1:La
            misc(k,i) = misc(k,i) + b(j,i)*diff(a(k),y(j));
        end
    end
end
L1a = misc;
fprintf("L1a is created.\n")

%% Ljb
    misc = sym(zeros(size(b)));
for i = 1:Lb    
    for k = 1:La
        for j = 1:La
            misc(k,i) = misc(k,i) + b(j,i)*diff(b(k,i),y(j));
        end
    end
end
L1b = misc;
fprintf("L1b is created.\n")

%% L1L1b
    misc = sym(zeros(size(b)));
for i = 1:Lb    
    for k = 1:La
        for j = 1:La
            misc(k,i) = misc(k,i) + b(j,i)*diff(L1b(k,i),y(j));
        end
    end
end
L1L1b = misc;
fprintf("L1L1b is created.\n")

end

