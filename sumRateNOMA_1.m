function SR_NOMA=sumRateNOMA_1(Pt,Pr,h1,h2,h3,g1,g2,NT,p_1_ini,p_2_ini,p_c_ini,theta_ini,ind_relay,tolerance,num_iter)
% user 2 as the common message,p_c=p_2,p_2=0
%% introduced parameters
T_1=abs(h1*p_1_ini)^2+1;
T_2=abs(h2*p_1_ini)^2+1;
T_3=abs(g1*p_1_ini)^2+1;
T_4=abs(g1*p_1_ini)^2;
% T_5=0;
T_c_1=abs(h1*p_2_ini)^2;
T_c_2=abs(h2*p_2_ini)^2;
T_c_3=abs(g1*p_2_ini)^2;
I_1=T_1-abs(h1*p_1_ini)^2;
% I_2=T_2;
I_3=T_3-abs(g1*p_1_ini)^2+T_c_3;
% I_4=T_3+T_c_3;

% \rho_c=[rho_c_1_ini,rho_c_2_ini,rho_c_e_ini]
rho_c_1_ini=T_c_1/T_1;
rho_c_2_ini=T_c_2/T_2;
rho_c_e_ini=T_c_3/T_3;

% \rho_p=[rho_p_1_ini,rho_p_2_ini,rho_1_e_ini,rho_2_e_ini]
rho_p_1_ini=T_1/I_1-1;
% rho_p_2_ini=T_2/I_2-1;
rho_1_e_ini=T_4/I_3;
% rho_2_e_ini=T_5/I_4;

beta_c_1_ini=log2(1+rho_c_1_ini);
beta_c_2_ini=log2(1+rho_c_2_ini);
beta_c_e_ini=log2(1+rho_c_e_ini);

beta_p_1_ini=log2(1+rho_p_1_ini);
% beta_p_2_ini=log2(1+rho_p_2_ini);
beta_1_e_ini=log2(1+rho_1_e_ini);
% beta_2_e_ini=log2(1+rho_2_e_ini);

p_1=sdpvar(NT,1,'full','complex');
p_2=sdpvar(NT,1,'full','complex');
% p_c=sdpvar(NT,1,'full','complex');
theta=sdpvar;
alpha_p_1=sdpvar;
% alpha_p_2=sdpvar;
alpha_1_e=sdpvar;
% alpha_2_e=sdpvar;   
alpha_c_1=sdpvar;
alpha_c_2=sdpvar;
alpha_c_e=sdpvar;

beta_p_1=sdpvar;
% beta_p_2=sdpvar;
beta_1_e=sdpvar;
% beta_2_e=sdpvar;   
beta_c_1=sdpvar;
beta_c_2=sdpvar;
beta_c_e=sdpvar;

rho_p_1=sdpvar;
% rho_p_2=sdpvar;
rho_1_e=sdpvar;
% rho_2_e=sdpvar;   
rho_c_1=sdpvar;
rho_c_2=sdpvar;
rho_c_e=sdpvar;

Rc1=theta*log2(1+(h1*p_2)'*(h1*p_2)/((h1*p_1)'*(h1*p_1)+1));
Rc2=theta*log2(1+(h2*p_2)'*(h2*p_2)/((h2*p_1)'*(h2*p_1)+1))+(1-theta)*log2(1+Pr*h3'*h3);
Cce=theta*log2(1+(g1*p_2)'*(g1*p_2)/((g1*p_1)'*(g1*p_1)+1))+(1-theta)*log2(1+Pr*g2'*g2);
Rp1=theta*log2(1+(h1*p_1)'*(h1*p_1));
C1e=theta*log2(1+(g1*p_1)'*(g1*p_1)/((g1*p_2)'*(g1*p_2)+1));
Rc_Cce=min(Rc1,Rc2)-Cce;
Rp1_C1e=Rp1-C1e;

%  ops = sdpsettings('solver','fmincon','fmincon.Algorithm','sqp','fmincon.MaxIter',Inf,'fmincon.MaxFunEvals',Inf,'verbose',0);
% ops = sdpsettings('solver','fmincon','fmincon.MaxIter',Inf,'fmincon.MaxFunEvals',Inf,'verbose',1);
ops = sdpsettings('verbose',0);
% Objective=min(alpha_c_1,alpha_c_2)+alpha_p_1-alpha_1_e+alpha_p_2-alpha_2_e-alpha_c_e;  
Objective=min(alpha_c_1,alpha_c_2)+alpha_p_1-alpha_1_e-alpha_c_e;  

loop=1;
t_past=0;
count=0;

while (loop)
      
    %constraints    
    constraints(1)=-theta;
    constraints(2)=theta-1;
    
    constraints(3)=-alpha_p_1;
%     constraints(4)=-alpha_p_2;
    constraints(5)=-alpha_1_e;
%     constraints(6)=-alpha_2_e;
    constraints(7)=-alpha_c_1;
    constraints(8)=-alpha_c_2;
    constraints(9)=-alpha_c_e;
    
    constraints(10)=-beta_p_1;
%     constraints(11)=-beta_p_2;
    constraints(12)=-beta_1_e;
%     constraints(13)=-beta_2_e;
    constraints(14)=-beta_c_1;
    constraints(15)=-beta_c_2;
    constraints(16)=-beta_c_e; 
    
    constraints(17)=-rho_p_1+1e-6;
%     constraints(18)=-rho_p_2;
    constraints(19)=-rho_1_e+1e-6;
%     constraints(20)=-rho_2_e;
    constraints(21)=-rho_c_1+1e-6;
    constraints(22)=-rho_c_2+1e-6;
    constraints(23)=-rho_c_e+1e-6;

%     constraints(24)=p_1'*p_1+p_2'*p_2+p_c'*p_c-Pt; % Quadratic scalar (real, 24 variables)
    constraints(24)=p_1'*p_1+p_2'*p_2-Pt; % Quadratic scalar (real, 24 variables)
    
    constraints(25)=-0.5*(theta_ini+beta_p_1_ini)*(theta+beta_p_1)+0.25*(theta_ini+beta_p_1_ini)^2+0.25*(theta-beta_p_1)^2+alpha_p_1; % Quadratic scalar (real, 3 variables)
%     constraints(26)=-0.5*(theta_ini+beta_p_2_ini)*(theta+beta_p_2)+0.25*(theta_ini+beta_p_2_ini)^2+0.25*(theta-beta_p_2)^2+alpha_p_2;

%     if ind_relay==1     
        constraints(27)=-0.5*(theta_ini+beta_c_1_ini)*(theta+beta_c_1)+0.25*(theta_ini+beta_c_1_ini)^2+0.25*(theta-beta_c_1)^2+alpha_c_1;
        constraints(28)=-0.5*(theta_ini+beta_c_2_ini)*(theta+beta_c_2)+0.25*(theta_ini+beta_c_2_ini)^2+0.25*(theta-beta_c_2)^2-(1-theta)*log2(1+Pr*(h3'*h3))+alpha_c_2;  
%      else
%         constraints(27)=-0.5*(theta_ini+beta_c_1_ini)*(theta+beta_c_1)+0.25*(theta_ini+beta_c_1_ini)^2+0.25*(theta-beta_c_1)^2-(1-theta)*log2(1+Pr*(h3'*h3))+alpha_c_1;
%         constraints(28)=-0.5*(theta_ini+beta_c_2_ini)*(theta+beta_c_2)+0.25*(theta_ini+beta_c_2_ini)^2+0.25*(theta-beta_c_2)^2+alpha_c_2;
%      end
    constraints(29)=-0.5*(theta_ini-beta_1_e_ini)*(theta-beta_1_e)+0.25*(theta_ini-beta_1_e_ini)^2+0.25*(theta+beta_1_e)^2-alpha_1_e;
%     constraints(30)=-0.5*(theta_ini-beta_2_e_ini)*(theta-beta_2_e)+0.25*(theta_ini-beta_2_e_ini)^2+0.25*(theta+beta_2_e)^2-alpha_2_e;
    constraints(31)=-0.5*(theta_ini-beta_c_e_ini)*(theta-beta_c_e)+0.25*(theta_ini-beta_c_e_ini)^2+0.25*(theta+beta_c_e)^2-alpha_c_e+(1-theta)*log2(1+Pr*(g2'*g2));

%     constraints(32)=alpha_c_e-alpha_c_1; % Linear scalar (real, 2 variables)
%     constraints(33)=alpha_c_e-alpha_c_2;
    
    % Quadratic scalar (real)
% %    constraints(34)=((h1*p_2)'*(h1*p_2))+1-2*real(p_1_ini'*h1'*h1*p_1)/rho_p_1_ini+((h1*p_1_ini)'*(h1*p_1_ini))*rho_p_1/(rho_p_1_ini^2);
    constraints(34)=1-2*real(p_1_ini'*h1'*h1*p_1)/rho_p_1_ini+((h1*p_1_ini)'*(h1*p_1_ini))*rho_p_1/(rho_p_1_ini^2);
%     constraints(35)=((h2*p_1)'*(h2*p_1))+1-2*real(p_2_ini'*h2'*h2*p_2)/rho_p_2_ini+((h2*p_2_ini)'*(h2*p_2_ini))*rho_p_2/(rho_p_2_ini^2);  
% %    constraints(36)=((h1*p_1)'*(h1*p_1))+((h1*p_2)'*(h1*p_2))+1-2*real(p_c_ini'*h1'*h1*p_c)/rho_c_1_ini+((h1*p_c_ini)'*(h1*p_c_ini))*rho_c_1/(rho_c_1_ini^2);
% %    constraints(37)=((h2*p_1)'*(h2*p_1))+((h2*p_2)'*(h2*p_2))+1-2*real(p_c_ini'*h2'*h2*p_c)/rho_c_2_ini+((h2*p_c_ini)'*(h2*p_c_ini))*rho_c_2/(rho_c_2_ini^2);

    constraints(36)=((h1*p_1)'*(h1*p_1))+1-2*real(p_2_ini'*h1'*h1*p_2)/rho_c_1_ini+((h1*p_2_ini)'*(h1*p_2_ini))*rho_c_1/(rho_c_1_ini^2);
    constraints(37)=((h2*p_1)'*(h2*p_1))+1-2*real(p_2_ini'*h2'*h2*p_2)/rho_c_2_ini+((h2*p_2_ini)'*(h2*p_2_ini))*rho_c_2/(rho_c_2_ini^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% rewrite begin %%%%%%%%%%%%%%%%
    % Signomial scalar(real)
%     constraints(38)=-rho_1_e*(2*real(p_c_ini'*g1'*g1*p_c+p_2_ini'*g1'*g1*p_2)-((g1*p_c_ini)'*(g1*p_c_ini))-((g1*p_2_ini)'*(g1*p_2_ini))+1)+((g1*p_1)'*(g1*p_1));%/rho_1_e;  
%     constraints(39)=-rho_2_e*(2*real(p_c_ini'*g1'*g1*p_c+p_1_ini'*g1'*g1*p_1)-((g1*p_c_ini)'*(g1*p_c_ini))-((g1*p_1_ini)'*(g1*p_1_ini))+1)+((g1*p_2)'*(g1*p_2));%/rho_2_e; 
%     constraints(40)=-rho_c_e*(2*real(p_1_ini'*g1'*g1*p_1+p_2_ini'*g1'*g1*p_2)-((g1*p_1_ini)'*(g1*p_1_ini))-((g1*p_2_ini)'*(g1*p_2_ini))+1)+((g1*p_c)'*(g1*p_c));%/rho_c_e;  
%     constraints(38)=-rho_1_e*((g1*p_c)*(g1*p_c)'+(g1*p_2)*(g1*p_2)'+1)+(g1*p_1)*(g1*p_1)';
%     constraints(39)=-rho_2_e*((g1*p_c)*(g1*p_c)'+(g1*p_1)*(g1*p_1)'+1)+(g1*p_2)*(g1*p_2)';
%     constraints(40)=-rho_c_e*((g1*p_1)*(g1*p_1)'+(g1*p_2)*(g1*p_2)'+1)+(g1*p_c)*(g1*p_c)';

% %    constraints(38)=-(((g1*p_c_ini)'*(g1*p_c_ini)+(g1*p_2_ini)'*(g1*p_2_ini))*(rho_1_e-2*rho_1_e_ini)+rho_1_e+2*rho_1_e_ini*real(p_c_ini'*g1'*g1*p_c+p_2_ini'*g1'*g1*p_2))+(g1*p_1)'*(g1*p_1);
%     constraints(39)=-(((g1*p_c_ini)'*(g1*p_c_ini)+(g1*p_1_ini)'*(g1*p_1_ini))*(rho_2_e-2*rho_2_e_ini)+rho_2_e+2*rho_2_e_ini*real(p_c_ini'*g1'*g1*p_c+p_1_ini'*g1'*g1*p_1))+(g1*p_2)'*(g1*p_2);
% %    constraints(40)=-(((g1*p_1_ini)'*(g1*p_1_ini)+(g1*p_2_ini)'*(g1*p_2_ini))*(rho_c_e-2*rho_c_e_ini)+rho_c_e+2*rho_c_e_ini*real(p_1_ini'*g1'*g1*p_1+p_2_ini'*g1'*g1*p_2))+(g1*p_c)'*(g1*p_c);
    constraints(38)=-(((g1*p_2_ini)'*(g1*p_2_ini))*(rho_1_e-2*rho_1_e_ini)+rho_1_e+2*rho_1_e_ini*real(p_2_ini'*g1'*g1*p_2))+(g1*p_1)'*(g1*p_1);
    constraints(40)=-(((g1*p_1_ini)'*(g1*p_1_ini))*(rho_c_e-2*rho_c_e_ini)+rho_c_e+2*rho_c_e_ini*real(p_1_ini'*g1'*g1*p_1))+(g1*p_2)'*(g1*p_2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% rewrite end %%%%%%%%%%%%%%%%%%%%%%%
    % Linear scalar (real)   
    constraints(41)=1+rho_1_e-(1+log(2)*(beta_1_e-beta_1_e_ini))*exp(beta_1_e_ini*log(2));
%     constraints(42)=1+rho_2_e-(1+log(2)*(beta_2_e-beta_2_e_ini))*exp(beta_2_e_ini*log(2));
    constraints(43)=1+rho_c_e-(1+log(2)*(beta_c_e-beta_c_e_ini))*exp(beta_c_e_ini*log(2));   
    
    % Nonlinear scalar (real)
    constraints(44)=exp(beta_p_1*log(2))-(1+rho_p_1); 
%     constraints(45)=exp(beta_p_2*log(2))-(1+rho_p_2);
    constraints(46)=exp(beta_c_1*log(2))-(1+rho_c_1);
    constraints(47)=exp(beta_c_2*log(2))-(1+rho_c_2);

    Constraints=[constraints(1:3)<=0,...
                constraints(5)<=0,...
                constraints(7:10)<=0,...
                constraints(12)<=0,...
                constraints(14:17)<=0,...
                constraints(19)<=0,...
                constraints(21:25)<=0,...
                constraints(27:29)<=0,...
                constraints(31)<=0,...
                constraints(34)<=0,...
                constraints(36:37)<=0,...
                constraints(38)<=0,...
                constraints(40:41)<=0,...
                constraints(43:44)<=0,...
                constraints(46:47)<=0];
    
    P = optimize(Constraints,-Objective,ops);

    if abs(value(Objective)-t_past)<=tolerance
        loop=0;
    else
        t_past=value(Objective);
        count=count+1;
        p_1_ini=value(p_1);
%         p_2_ini=value(p_2);
        p_2_ini=value(p_2);
        theta_ini=value(theta);
        beta_p_1_ini=value(beta_p_1);
%         beta_p_2_ini=value(beta_p_2);
        beta_1_e_ini=value(beta_1_e);
%         beta_2_e_ini=value(beta_2_e);
        beta_c_1_ini=value(beta_c_1);
        beta_c_2_ini=value(beta_c_2);
        beta_c_e_ini=value(beta_c_e);

        rho_p_1_ini=value(rho_p_1);
%         rho_p_2_ini=value(rho_p_2);
        rho_1_e_ini=value(rho_1_e);
%         rho_2_e_ini=value(rho_2_e);
        rho_c_1_ini=value(rho_c_1);
        rho_c_2_ini=value(rho_c_2);
        rho_c_e_ini=value(rho_c_e);
    end
    
    if count>=num_iter
        break;
    end
   
end

%% show the outputs
value_Rc_Cce=value(Rc_Cce);
value_Rp1_C1e=value(Rp1_C1e);

if value_Rc_Cce<=0 && value_Rp1_C1e<=0
    SR_NOMA=0;
end
if value_Rc_Cce>=0 && value_Rp1_C1e<=0
    SR_NOMA=value_Rc_Cce;
end
if value_Rc_Cce<=0 && value_Rp1_C1e>=0
    SR_NOMA=value_Rp1_C1e;
end
if value_Rc_Cce>=0 && value_Rp1_C1e>=0
    SR_NOMA=value_Rc_Cce+value_Rp1_C1e;
end
theta_opt=value(theta);

fprintf('NOMA_1 : theta= %1.0f,rate= %1.3E \n',[theta_opt SR_NOMA]);
end