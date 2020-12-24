function [SR_RS_4]=sumRateRS_4(Pt,Pr,h1,h2,h3,g1,g2,NT,p_1_ini,p_2_ini,p_c_ini,theta_ini,ind_relay,tolerance,num_iter,epsilon_1,epsilon_2)
%% introduced parameters
% 计算初始值的时候 先用g1,g2的估计值来代替g1,g2
T_1=abs(h1*p_1_ini)^2+abs(h1*p_2_ini)^2+1;
T_2=abs(h2*p_1_ini)^2+abs(h2*p_2_ini)^2+1;
T_3=abs(g1*p_1_ini)^2+abs(g1*p_2_ini)^2+1;
T_4=abs(g1*p_1_ini)^2;
T_5=abs(g1*p_2_ini)^2;
T_c_1=abs(h1*p_c_ini)^2;
T_c_2=abs(h2*p_c_ini)^2;
T_c_3=abs(g1*p_c_ini)^2;
I_1=T_1-abs(h1*p_1_ini)^2;
I_2=T_2-abs(h2*p_2_ini)^2;
I_3=T_3-abs(g1*p_1_ini)^2+T_c_3;
I_4=T_3-abs(g1*p_2_ini)^2+T_c_3;

% \rho_c=[rho_c_1_ini,rho_c_2_ini,rho_c_e_ini_1,rho_c_e_ini_2]
rho_c_1_ini=T_c_1/T_1;
rho_c_2_ini=T_c_2/T_2;
rho_c_e_ini_1=T_c_3/T_3;
rho_c_e_ini_2=Pr*abs(g2)^2;


% \rho_p=[rho_p_1_ini,rho_p_2_ini,rho_1_e_ini,rho_2_e_ini]
rho_p_1_ini=T_1/I_1-1;
rho_p_2_ini=T_2/I_2-1;
rho_1_e_ini=T_4/I_3;
rho_2_e_ini=T_5/I_4;

beta_c_1_ini=log2(1+rho_c_1_ini);
beta_c_2_ini=log2(1+rho_c_2_ini);
beta_c_e_ini_1=log2(1+rho_c_e_ini_1);
beta_c_e_ini_2=log2(1+rho_c_e_ini_2);

beta_p_1_ini=log2(1+rho_p_1_ini);
beta_p_2_ini=log2(1+rho_p_2_ini);
beta_1_e_ini=log2(1+rho_1_e_ini);
beta_2_e_ini=log2(1+rho_2_e_ini);

m_1_ini = abs(h1*p_c_ini)^2/rho_c_1_ini;
m_2_ini = abs(h2*p_c_ini)^2/rho_c_2_ini;
n_1_ini = abs(h1*p_1_ini)^2/rho_p_1_ini;
n_2_ini = abs(h2*p_2_ini)^2/rho_p_2_ini;

d_ini = abs(g1*p_c_ini)^2/rho_c_e_ini_1;
e_1_ini = abs(g1*p_1_ini)^2/rho_1_e_ini;
e_2_ini = abs(g1*p_2_ini)^2/rho_2_e_ini;

P_1=sdpvar(NT,NT,'full','complex');
P_2=sdpvar(NT,NT,'full','complex');
P_c=sdpvar(NT,NT,'full','complex');

t=sdpvar;
theta=sdpvar;
alpha_p_1=sdpvar;
alpha_p_2=sdpvar;
alpha_1_e=sdpvar;
alpha_2_e=sdpvar;   
alpha_c_1=sdpvar;
alpha_c_2=sdpvar;
alpha_c_e=sdpvar;

beta_p_1=sdpvar;
beta_p_2=sdpvar;
beta_1_e=sdpvar;
beta_2_e=sdpvar;   
beta_c_1=sdpvar;
beta_c_2=sdpvar;
beta_c_e_1=sdpvar;
beta_c_e_2=sdpvar;

rho_p_1=sdpvar;
rho_p_2=sdpvar;
rho_1_e=sdpvar;
rho_2_e=sdpvar;   
rho_c_1=sdpvar;
rho_c_2=sdpvar;
rho_c_e_1=sdpvar;
rho_c_e_2=sdpvar;

delta_1 = sdpvar;
delta_2 = sdpvar;
delta_3 = sdpvar;
delta_4 = sdpvar;
delta_5 = sdpvar;
delta_6 = sdpvar;
delta_7 = sdpvar;

m_1 = sdpvar;
m_2 = sdpvar;
n_1 = sdpvar;
n_2 = sdpvar;

d = sdpvar;
e_1 = sdpvar;
e_2 = sdpvar;



%  ops = sdpsettings('solver','fmincon','fmincon.Algorithm','sqp','fmincon.MaxIter',Inf,'fmincon.MaxFunEvals',Inf,'verbose',0);
% ops = sdpsettings('solver','fmincon','fmincon.MaxIter',Inf,'fmincon.MaxFunEvals',Inf,'verbose',1);
ops = sdpsettings('solver','mosek','verbose',0,'debug',0);

Objective = t;  

loop=1;
t_past=0;
count=0;

while (loop)
    A = P_1+P_2;
    B_1 = P_c+P_2;
    B_2 = P_c+P_1;
    
    theta_1_e = -0.5*(theta_ini-beta_1_e_ini)*(theta-beta_1_e)+0.25*(theta_ini-beta_1_e_ini)^2+0.25*(theta+beta_1_e)^2;
    theta_2_e = -0.5*(theta_ini-beta_2_e_ini)*(theta-beta_2_e)+0.25*(theta_ini-beta_2_e_ini)^2+0.25*(theta+beta_2_e)^2;
    theta_c_e_1 = -0.5*(theta_ini-beta_c_e_ini_1)*(theta-beta_c_e_1)+0.25*(theta_ini-beta_c_e_ini_1)^2+0.25*(theta+beta_c_e_1)^2;
    theta_n1_rho=-0.5*(n_1_ini-rho_p_1_ini)*(n_1-rho_p_1)+0.25*(n_1_ini-rho_p_1_ini)^2+0.25*(n_1+rho_p_1)^2;
    theta_n2_rho=-0.5*(n_2_ini-rho_p_2_ini)*(n_2-rho_p_2)+0.25*(n_2_ini-rho_p_2_ini)^2+0.25*(n_2+rho_p_2)^2;
    theta_m1_rho=-0.5*(m_1_ini-rho_c_1_ini)*(m_1-rho_c_1)+0.25*(m_1_ini-rho_c_1_ini)^2+0.25*(m_1+rho_c_1)^2;
    theta_m2_rho=-0.5*(m_2_ini-rho_c_2_ini)*(m_2-rho_c_2)+0.25*(m_2_ini-rho_c_2_ini)^2+0.25*(m_2+rho_c_2)^2;
    
    phi_p_1 = 0.5*(theta_ini+beta_p_1_ini)*(theta+beta_p_1)-0.25*(theta_ini+beta_p_1_ini)^2-0.25*(theta-beta_p_1)^2;
    phi_p_2 = 0.5*(theta_ini+beta_p_2_ini)*(theta+beta_p_2)-0.25*(theta_ini+beta_p_2_ini)^2-0.25*(theta-beta_p_2)^2;
    phi_c_1 = 0.5*(theta_ini+beta_c_1_ini)*(theta+beta_c_1)-0.25*(theta_ini+beta_c_1_ini)^2-0.25*(theta-beta_c_1)^2;
    phi_c_2 = 0.5*(theta_ini+beta_c_2_ini)*(theta+beta_c_2)-0.25*(theta_ini+beta_c_2_ini)^2-0.25*(theta-beta_c_2)^2;
    phi_c_e_2 = 0.5*(theta_ini+beta_c_e_ini_2)*(theta+beta_c_e_2)-0.25*(theta_ini+beta_c_e_ini_2)^2-0.25*(theta-beta_c_e_2)^2;
    phi_d_rho = 0.5*(d_ini+rho_c_e_ini_1)*(d+rho_c_e_1)-0.25*(d_ini+rho_c_e_ini_1)^2-0.25*(d-rho_c_e_1)^2;
    phi_e_rho_1 = 0.5*(e_1_ini+rho_1_e_ini)*(e_1+rho_1_e)-0.25*(e_1_ini+rho_1_e_ini)^2-0.25*(e_1-rho_1_e)^2;
    phi_e_rho_2 = 0.5*(e_2_ini+rho_2_e_ini)*(e_2+rho_2_e)-0.25*(e_2_ini+rho_2_e_ini)^2-0.25*(e_2-rho_2_e)^2;
      
    %constraints    
    constraints(1)=-theta;
    constraints(2)=theta-1;
    
    constraints(3)=-delta_1;
    constraints(4)=-delta_2;
    constraints(5)=-delta_3;
    constraints(6)=-delta_4;
    constraints(7)=-delta_5;
    constraints(8)=-delta_6;
    constraints(9)=-delta_7;
    
    constraints(10)=-m_1;
    constraints(11)=-m_2;
    constraints(12)=-n_1;
    constraints(13)=-n_2;
    
    constraints(14)=-d;
    constraints(15)=-e_1;
    constraints(16)=-e_2;
     
    constraints(17)=-alpha_p_1;
    constraints(18)=-alpha_p_2;
    constraints(19)=-alpha_1_e;
    constraints(20)=-alpha_2_e;
    constraints(21)=-alpha_c_1;
    constraints(22)=-alpha_c_2;
    constraints(23)=-alpha_c_e;
    
    constraints(24)=-beta_p_1;
    constraints(25)=-beta_p_2;
    constraints(26)=-beta_1_e;
    constraints(27)=-beta_2_e;
    constraints(28)=-beta_c_1;
    constraints(29)=-beta_c_2;
    constraints(30)=-beta_c_e_1; 
    constraints(31)=-beta_c_e_2; 
    
    constraints(32)=-rho_p_1+1e-4;
    constraints(33)=-rho_p_2+1e-4;
    constraints(34)=-rho_1_e+1e-4;
    constraints(35)=-rho_2_e+1e-4;
    constraints(36)=-rho_c_1+1e-4;
    constraints(37)=-rho_c_2+1e-4;
    constraints(38)=-rho_c_e_1+1e-4;
    constraints(39)=-rho_c_e_2+1e-4;

    constraints(40)=trace(P_1)+trace(P_2)+trace(P_c)-Pt; % Quadratic scalar 

    constraints(41)=-phi_p_1+alpha_p_1; % Quadratic scalar (real, 3 variables)
    constraints(42)=-phi_p_2+alpha_p_2;

    if ind_relay==1     
        constraints(43)=-phi_c_1+alpha_c_1;% Quadratic scalar 
        constraints(44)=-phi_c_2-(1-theta)*log2(1+Pr*(h3'*h3))+alpha_c_2;  
     else
        constraints(43)=-phi_c_1-(1-theta)*log2(1+Pr*(h3'*h3))+alpha_c_1;
        constraints(44)=-phi_c_2+alpha_c_2;
     end
    constraints(45)=theta_1_e-alpha_1_e;
    constraints(46)=theta_2_e-alpha_2_e;
    constraints(47)=theta_c_e_1-phi_c_e_2+beta_c_e_2-alpha_c_e;

    constraints(48)=alpha_c_e-alpha_c_1; % Linear scalar (real, 2 variables)
    constraints(49)=alpha_c_e-alpha_c_2;
    
    % Quadratic scalar (real)
    constraints(50)= h1*P_2*h1'+1-n_1;
    constraints(51)= theta_n1_rho-h1*P_1*h1';
    constraints(52)= h2*P_1*h2'+1-n_2;
    constraints(53)= theta_n2_rho-h2*P_2*h2';
    constraints(54)= h1*(P_1+P_2)*h1'+1-m_1;
    constraints(55)= theta_m1_rho-h1*P_c*h1';
    constraints(56)= h2*(P_1+P_2)*h2'+1-m_2;
    constraints(57)= theta_m2_rho-h2*P_c*h2';

    % Linear scalar (real)   
    constraints(58)=1+rho_1_e-(1+log(2)*(beta_1_e-beta_1_e_ini))*exp(beta_1_e_ini*log(2));
    constraints(59)=1+rho_2_e-(1+log(2)*(beta_2_e-beta_2_e_ini))*exp(beta_2_e_ini*log(2));
    constraints(60)=1+rho_c_e_1-(1+log(2)*(beta_c_e_1-beta_c_e_ini_1))*exp(beta_c_e_ini_1*log(2));  
    constraints(61)=1+rho_c_e_2-(1+log(2)*(beta_c_e_2-beta_c_e_ini_2))*exp(beta_c_e_ini_2*log(2)); 
    
%      Nonlinear scalar (real)
    constraints(62)=exp(beta_p_1*log(2))-(1+rho_p_1); 
    constraints(63)=exp(beta_p_2*log(2))-(1+rho_p_2);
    constraints(64)=exp(beta_c_1*log(2))-(1+rho_c_1);
    constraints(65)=exp(beta_c_2*log(2))-(1+rho_c_2);
    
    constraints(66)=t-(min(alpha_c_1,alpha_c_2)-alpha_c_e+alpha_p_1-alpha_1_e);
    constraints(67)=t-(min(alpha_c_1,alpha_c_2)-alpha_c_e+alpha_p_2-alpha_2_e);

    
    mat1=[P_c-delta_1*eye(NT),    P_c*g1'; 
          g1*P_c,                 g1*P_c*g1'-phi_d_rho+delta_1*epsilon_1^2];
      
    mat2=[A+delta_2*eye(NT),           A*g1'; 
          g1*A,                        g1*A*g1'+1-delta_2*epsilon_1^2-d];

    mat3=[Pr-delta_3,                  g2'*Pr; 
          g2*Pr,                       g2*g2'*Pr-rho_c_e_2+delta_3*epsilon_2^2];
      
    mat4=[P_1-delta_4*eye(NT),    P_1*g1'; 
          g1*P_1,                 g1*P_1*g1'-phi_e_rho_1+delta_4*epsilon_1^2];
      
    mat5=[P_2-delta_5*eye(NT),    P_2*g1'; 
            g1*P_2,               g1*P_2*g1'-phi_e_rho_2+delta_5*epsilon_1^2];  

    mat6=[B_1+delta_6*eye(NT),         B_1*g1'; 
          g1*B_1,                      g1*B_1*g1'+1-delta_6*epsilon_1^2-e_1];
      
    mat7=[B_2+delta_7*eye(NT),         B_2*g1'; 
          g1*B_2,                      g1*B_2*g1'+1-delta_7*epsilon_1^2-e_2];

    Constraints=[constraints(1:67)<=0,...
                P_1 >= 0,...
                P_2 >= 0,...
                P_c >= 0,...
                mat1<=0,...
                mat2>=0,...
                mat3<=0,...
                mat4<=0,...
                mat5<=0,...
                mat6>=0,...
                mat7>=0];
%     Constraints=[constraints(1:35)<=0,mat1<=0,mat2>=0,mat3<=0,mat4<=0,mat5<=0,mat6>=0,mat7>=0];
    
    P = optimize(Constraints,-Objective,ops);
%     P = optimize(Constraints,[],ops);    

    if abs(value(Objective)-t_past)<=tolerance
        loop=0;
    else
        t_past=value(Objective);
        count=count+1;

        theta_ini=value(theta);
        beta_p_1_ini=value(beta_p_1);
        beta_p_2_ini=value(beta_p_2);
        beta_1_e_ini=value(beta_1_e);
        beta_2_e_ini=value(beta_2_e);
        beta_c_1_ini=value(beta_c_1);
        beta_c_2_ini=value(beta_c_2);
        beta_c_e_ini_1=value(beta_c_e_1);
        beta_c_e_ini_2=value(beta_c_e_2);

        rho_p_1_ini=value(rho_p_1);
        rho_p_2_ini=value(rho_p_2);
        rho_1_e_ini=value(rho_1_e);
        rho_2_e_ini=value(rho_2_e);
        rho_c_1_ini=value(rho_c_1);
        rho_c_2_ini=value(rho_c_2);
        rho_c_e_ini_1=value(rho_c_e_1);
        rho_c_e_ini_2=value(rho_c_e_2);
        
        m_1_ini = value(m_1);
        m_2_ini = value(m_2);
        n_1_ini = value(n_1);
        n_2_ini = value(n_2);
        
        d_ini = value(d);
        e_1_ini = value(e_1);
        e_2_ini = value(e_2);
    end
%     
%     if rho_p_1_ini<=1e-4
%         break;
%     end
%     if rho_p_2_ini<=1e-4
%         break;
%     end
    
    if count>=num_iter
        break;
    end
end

%% show the outputs
% value_Rc_Cce=value(Rc_Cce);
% value_Rp1_C1e=value(Rp1_C1e);
% value_Rp2_C2e=value(Rp2_C2e);

% if value_Rp1_C1e<=0 && value_Rp2_C2e<=0
%     SR_RS_4=value_Rc_Cce;
% end
% if value_Rp1_C1e>=0 && value_Rp2_C2e<=0
%     SR_RS_4=value_Rc_Cce+value_Rp1_C1e;
% end
% if value_Rp1_C1e<=0 && value_Rp2_C2e>=0
%     SR_RS_4=value_Rc_Cce+value_Rp2_C2e;
% end
% if value_Rp1_C1e>=0 && value_Rp2_C2e>=0
%     SR_RS_4=value_Rc_Cce+value_Rp1_C1e+value_Rp2_C2e;
% end
SR_RS_4=value(t);
ratio_comm=value(trace(P_c))/Pt;
% ratio_comm=value(P_c'*P_c)/Pt;
theta_opt=value(theta);

fprintf('RS_4 : theta_4= %1.3f,rate_4= %1.3E,ratio_4=%1.3f \n',[theta_opt SR_RS_4 ratio_comm]);
end