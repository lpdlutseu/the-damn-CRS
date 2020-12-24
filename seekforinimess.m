clc;clear all; clf
%% parameter settings
db=25;
Pt=10^(db/10);
Pr=Pt;

bias_1=1;bias_2=1;bias_3=1;bias_e_1=5;bias_e_2=1;
epsilon_1 = 0.1; epsilon_2 = 0.1;

NT=4; N_user=2;  A_U=1;  tolerance = 10^-3;
randn('seed',2*N_user)

h1=sqrt(bias_1)/sqrt(2)*(randn(1,NT)+1i*randn(1,NT));
h2=sqrt(bias_2)/sqrt(2)*(randn(1,NT)+1i*randn(1,NT));
g1=sqrt(bias_e_1)/sqrt(2)*(randn(1,NT)+1i*randn(1,NT));

h3=norm(sqrt(bias_3)/sqrt(2)*(randn(1,NT)+1i*randn(1,NT)));
g2=norm(sqrt(bias_e_2)/sqrt(2)*(randn(1,NT)+1i*randn(1,NT)));

if norm(h1)>=norm(h2)
    ind_relay=1;
else
    ind_relay=2;
end
%% initial points way-1
p_1_ini=sqrt(Pt/N_user)*h1'/norm(h1)
p_2_ini=sqrt(Pt/N_user)*h2'/norm(h2)
p_c_ini=sqrt(Pt/N_user)*(h1+h2)'/norm(h1+h2)
theta_ini=0.5
%% initial points way-2
%     N_legi=N_user-1;
%     P_common=Pt*0.8;
%     P_private_k=(Pt-P_common)/N_legi;
% 
%     H=null(g1);
%     hat_p_c=H(:,1);
%     p_c_ini=hat_p_c/norm(hat_p_c)*sqrt(P_common);
%     p_1_ini=h1'/norm(h1)*sqrt(P_private_k);
%     p_2_ini=h2'/norm(h2)*sqrt(P_private_k);
%     theta_ini=0.5;
%% initial points way-3
%     N_legi=N_user-1;
%     P_common=Pt*0.8;
%     P_private_k=(Pt-P_common)/N_legi;
%     
%     H_tot=[h1' h2'];
%     [U2,~,~]=svd(H_tot);
%     hat_p_c=U2(:,1);
%     p_c_ini=hat_p_c*sqrt(P_common);
%     p_1_ini=h1'/norm(h1)*sqrt(P_private_k);
%     p_2_ini=h2'/norm(h2)*sqrt(P_private_k);
%     theta_ini=0.5;
%% optimization
p_1=sdpvar(NT,1,'full','complex');
p_2=sdpvar(NT,1,'full','complex');
p_c=sdpvar(NT,1,'full','complex');
theta=sdpvar;
assign(p_1,p_1_ini);
assign(p_2,p_2_ini);
assign(p_c,p_c_ini);
assign(theta,theta_ini);
ops=sdpsettings('verbose',0,'usex0',1);

%     if ind_relay==1
%         Objective=min(theta*log2(1+((h1*p_c)'*(h1*p_c))/(((h1*p_1)'*(h1*p_1))+((h1*p_2)'*(h1*p_2))+1)),...
%         theta*log2(1+((h2*p_c)'*(h2*p_c))/(((h2*p_1)'*(h2*p_1))+((h2*p_2)'*(h2*p_2))+1))+(1-theta)*log2(1+Pr*(h3'*h3)))-...
%         (theta*log2(1+((g1*p_c)'*(g1*p_c))/(((g1*p_1)'*(g1*p_1))+((g1*p_2)'*(g1*p_2))+1))+(1-theta)*log2(1+Pr*(g2'*g2)));
%     else
%         Objective=min(theta*log2(1+((h1*p_c)'*(h1*p_c))/(((h1*p_1)'*(h1*p_1))+((h1*p_2)'*(h1*p_2))+1))+(1-theta)*log2(1+Pr*(h3'*h3)),...
%         theta*log2(1+((h2*p_c)'*(h2*p_c))/(((h2*p_1)'*(h2*p_1))+((h2*p_2)'*(h2*p_2))+1)))-...
%         (theta*log2(1+((g1*p_c)'*(g1*p_c))/(((g1*p_1)'*(g1*p_1))+((g1*p_2)'*(g1*p_2))+1))+(1-theta)*log2(1+Pr*(g2'*g2)));
%     end 
%     Constraints=[p_1'*p_1+p_2'*p_2+p_c'*p_c-Pt<=0,...
%                      0 <= theta <= 1];

if ind_relay==1
    Constraints=[p_1'*p_1+p_2'*p_2+p_c'*p_c-Pt<=0,...
                 0 <= theta <= 1,...
                min(theta*log2(1+((h1*p_c)'*(h1*p_c))/(((h1*p_1)'*(h1*p_1))+((h1*p_2)'*(h1*p_2))+1)),...
                theta*log2(1+((h2*p_c)'*(h2*p_c))/(((h2*p_1)'*(h2*p_1))+((h2*p_2)'*(h2*p_2))+1))+(1-theta)*log2(1+Pr*(h3'*h3)))-...
                (theta*log2(1+((g1*p_c)'*(g1*p_c))/(((g1*p_1)'*(g1*p_1))+((g1*p_2)'*(g1*p_2))+1))+(1-theta)*log2(1+Pr*(g2'*g2)))>=0];
%         Constraints=[p_1'*p_1+p_2'*p_2+p_c'*p_c-Pt<=0,...
%                      0 <= theta <= 1,...
%                     min(theta*log2(1+((h1*p_c)'*(h1*p_c))/(((h1*p_1)'*(h1*p_1))+((h1*p_2)'*(h1*p_2))+1)),...
%                     theta*log2(1+((h2*p_c)'*(h2*p_c))/(((h2*p_1)'*(h2*p_1))+((h2*p_2)'*(h2*p_2))+1))+(1-theta)*log2(1+Pr*(h3'*h3)))-...
%                     (theta*log2(1+((g1*p_c)'*(g1*p_c))/(((g1*p_1)'*(g1*p_1))+((g1*p_2)'*(g1*p_2))+1))+(1-theta)*log2(1+Pr*(g2'*g2)))>=0,...
%                     min(theta*log2(1+((h1*p_1)'*(h1*p_1))/(((h1*p_2)'*(h1*p_2))+1)),...
%                     theta*log2(1+((h2*p_1)'*(h2*p_1))/(((h2*p_2)'*(h2*p_2))+1))+(1-theta)*log2(1+Pr*(h3'*h3)))-...
%                     (theta*log2(1+((g1*p_1)'*(g1*p_1))/(((g1*p_2)'*(g1*p_2))+1))+(1-theta)*log2(1+Pr*(g2'*g2)))>=0];
else
        Constraints=[p_1'*p_1+p_2'*p_2+p_c'*p_c-Pt<=0,...
                     0 <= theta <= 1,...
                     min(theta*log2(1+((h1*p_c)'*(h1*p_c))/(((h1*p_1)'*(h1*p_1))+((h1*p_2)'*(h1*p_2))+1))+(1-theta)*log2(1+Pr*(h3'*h3)),...
                    theta*log2(1+((h2*p_c)'*(h2*p_c))/(((h2*p_1)'*(h2*p_1))+((h2*p_2)'*(h2*p_2))+1)))-...
                    (theta*log2(1+((g1*p_c)'*(g1*p_c))/(((g1*p_1)'*(g1*p_1))+((g1*p_2)'*(g1*p_2))+1))+(1-theta)*log2(1+Pr*(g2'*g2)))>=0];
%             Constraints=[p_1'*p_1+p_2'*p_2+p_c'*p_c-Pt<=0,...
%                          0 <= theta <= 1,...
%                          min(theta*log2(1+((h1*p_c)'*(h1*p_c))/(((h1*p_1)'*(h1*p_1))+((h1*p_2)'*(h1*p_2))+1))+(1-theta)*log2(1+Pr*(h3'*h3)),...
%                         theta*log2(1+((h2*p_c)'*(h2*p_c))/(((h2*p_1)'*(h2*p_1))+((h2*p_2)'*(h2*p_2))+1)))-...
%                         (theta*log2(1+((g1*p_c)'*(g1*p_c))/(((g1*p_1)'*(g1*p_1))+((g1*p_2)'*(g1*p_2))+1))+(1-theta)*log2(1+Pr*(g2'*g2)))>=0,...
%                         min(theta*log2(1+((h1*p_2)'*(h1*p_2))/(((h1*p_1)'*(h1*p_1))+1))+(1-theta)*log2(1+Pr*(h3'*h3)),...
%                         theta*log2(1+((h2*p_2)'*(h2*p_2))/(((h2*p_1)'*(h2*p_1))+1)))-...
%                         (theta*log2(1+((g1*p_2)'*(g1*p_2))/(((g1*p_1)'*(g1*p_1))+1))+(1-theta)*log2(1+Pr*(g2'*g2)))>=0];
end    

%     P=optimize(Constraints,-Objective,ops);
P=optimize(Constraints,[],ops);

p_1=value(p_1)
p_2=value(p_2)
p_c=value(p_c)
theta=value(theta)
ratio_comm=value(p_c'*p_c)/Pt;
orth_g1pc=value(norm(g1*p_c));
fprintf('ini : theta= %1.0f,ratio=%1.3f,orth_g1pc=%1.3f \n',[theta ratio_comm orth_g1pc]);
