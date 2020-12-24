function [p_1_ini,p_2_ini,p_c_ini,theta_ini]=seek(Pt,Pr,h1,h2,h3,g1,g2,NT,N_user)

theta=0.5;
p_1=sqrt(Pt/N_user)*h1'/norm(h1);
p_2=sqrt(Pt/N_user)*h2'/norm(h2);
p_c=sqrt(Pt/N_user)*(h1+h2)'/norm(h1+h2);

% N_legi=N_user-1;
% H_tot=[h1' h2'];
% P_common=Pt*0.8;
% P_private_k=(Pt-P_common)/N_legi;
% 
% [U2,~,~]=svd(H_tot);
% hat_p_c=U2(:,1);
% 
% p_1=h1'/norm(h1)*sqrt(P_private_k);
% p_2=h2'/norm(h2)*sqrt(P_private_k);
% p_c=hat_p_c*sqrt(P_common);
% theta=0.8;

x0p12ct=[p_1;p_2;p_c;theta];
save('data.mat');
Amar=[zeros(1,3*NT),-1;zeros(1,3*NT),1];
buneq=[0;1];
x = fmincon(@objfun,x0p12ct,Amar,buneq,[],[],[],[],@nonlcon,optimoptions('fmincon','MaxFunEvals',Inf));

objfun(x)
p_1_ini=x(1:NT);
p_2_ini=x(NT+1:2*NT);
p_c_ini=x(NT*2+1:3*NT);
theta_ini=x(end);
ratio_comm=(p_c_ini'*p_c_ini)/Pt;
orth_g1pc=norm(g1*p_c_ini);
fprintf('ini : theta= %1.0f,ratio=%1.3f,orth_g1pc=%1.3f \n',[theta_ini ratio_comm orth_g1pc]);
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

function [c,ceq]=nonlcon(x)
load data

p_1=x(1:NT);
p_2=x(NT+1:2*NT);
p_c=x(NT*2+1:3*NT);
theta=x(end);
c=(p_1'*p_1)+(p_2'*p_2)+(p_c'*p_c)-Pt;
ceq=[];
end



end
