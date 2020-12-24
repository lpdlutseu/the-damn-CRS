function [c,ceq]=nonlcon(x)
load data

p_1=x(1:NT);
p_2=x(NT+1:2*NT);
p_c=x(NT*2+1:3*NT);
theta=x(end);
c=(p_1'*p_1)+(p_2'*p_2)+(p_c'*p_c)-Pt;
ceq=[];
end