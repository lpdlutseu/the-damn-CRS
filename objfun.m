function obj=objfun(x)
load data

p_1=x(1:NT);
p_2=x(NT+1:2*NT);
p_c=x(NT*2+1:3*NT);
theta=x(end);

obj=-min(theta*log2(1+(abs(h1*p_c)^2)/((abs(h1*p_1)^2)+(abs(h1*p_2)^2)+1)),...
    theta*log2(1+(abs(h2*p_c)^2)/((abs(h2*p_1)^2)+(abs(h2*p_2)^2)+1))+(1-theta)*log2(1+Pr*abs(h3)^2))+...
    (theta*log2(1+(abs(g1*p_c)^2)/((abs(g1*p_1)^2)+(abs(g1*p_2)^2)+1))+(1-theta)*log2(1+Pr*abs(g2)^2));
end