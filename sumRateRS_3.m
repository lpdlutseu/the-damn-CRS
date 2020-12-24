function SR_RS_3=sumRateRS_3(Pt,Pr,h1,h2,h3,g1,g2,NT,p_1_ini,p_2_ini,p_c_ini,theta_ini,ind_relay,tolerance,num_iter)
%% introduced parameters
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

% \rho_c=[rho_c_1_ini,rho_c_2_ini,rho_c_e_ini]
rho_c_1_ini=T_c_1/T_1;
rho_c_2_ini=T_c_2/T_2;
rho_c_e_ini=T_c_3/T_3;

% \rho_p=[rho_p_1_ini,rho_p_2_ini,rho_1_e_ini,rho_2_e_ini]
rho_p_1_ini=T_1/I_1-1;
rho_p_2_ini=T_2/I_2-1;
rho_1_e_ini=T_4/I_3;
rho_2_e_ini=T_5/I_4;

beta_c_1_ini=log2(1+rho_c_1_ini);
beta_c_2_ini=log2(1+rho_c_2_ini);
beta_c_e_ini=log2(1+rho_c_e_ini);

beta_p_1_ini=log2(1+rho_p_1_ini);
beta_p_2_ini=log2(1+rho_p_2_ini);
beta_1_e_ini=log2(1+rho_1_e_ini);
beta_2_e_ini=log2(1+rho_2_e_ini);

p_1=sdpvar(NT,1,'full','complex');
p_2=sdpvar(NT,1,'full','complex');
p_c=sdpvar(NT,1,'full','complex');
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
beta_c_e=sdpvar;

rho_p_1=sdpvar;
rho_p_2=sdpvar;
rho_1_e=sdpvar;
rho_2_e=sdpvar;   
rho_c_1=sdpvar;
rho_c_2=sdpvar;
rho_c_e=sdpvar;

% ops = sdpsettings('solver','fmincon','fmincon.Algorithm','sqp','fmincon.MaxIter',Inf,'fmincon.MaxFunEvals',Inf,'verbose',1);
% ops = sdpsettings('solver','fmincon','fmincon.MaxIter',Inf,'fmincon.MaxFunEvals',Inf,'verbose',1);
ops = sdpsettings('verbose',0);
Objective=min(alpha_c_1,alpha_c_2)-alpha_c_e+alpha_p_2-alpha_2_e;

loop=1;
t_past=0;
count=0;


while (loop)
      
    %constraints    
    constraints(1)=-theta;
    constraints(2)=theta-1;

    constraints(3)=-alpha_c_1;
    constraints(4)=-alpha_c_2;
    constraints(5)=-alpha_c_e;

    constraints(6)=-beta_c_1;
    constraints(7)=-beta_c_2;
    constraints(8)=-beta_c_e;

    constraints(9)=-rho_c_1;
    constraints(10)=-rho_c_2;
    constraints(11)=-rho_c_e;

    constraints(12)=p_1'*p_1+p_2'*p_2+p_c'*p_c-Pt; % Quadratic scalar (real, 24 variables)  

    if ind_relay==1     
        constraints(13)=-0.5*(theta_ini+beta_c_1_ini)*(theta+beta_c_1)+0.25*(theta_ini+beta_c_1_ini)^2+0.25*(theta-beta_c_1)^2+alpha_c_1;
        constraints(14)=-0.5*(theta_ini+beta_c_2_ini)*(theta+beta_c_2)+0.25*(theta_ini+beta_c_2_ini)^2+0.25*(theta-beta_c_2)^2-(1-theta)*log2(1+Pr*(h3'*h3))+alpha_c_2;  
     else
        constraints(13)=-0.5*(theta_ini+beta_c_1_ini)*(theta+beta_c_1)+0.25*(theta_ini+beta_c_1_ini)^2+0.25*(theta-beta_c_1)^2-(1-theta)*log2(1+Pr*(h3'*h3))+alpha_c_1;
        constraints(14)=-0.5*(theta_ini+beta_c_2_ini)*(theta+beta_c_2)+0.25*(theta_ini+beta_c_2_ini)^2+0.25*(theta-beta_c_2)^2+alpha_c_2;
    end

    constraints(15)=-0.5*(theta_ini-beta_c_e_ini)*(theta-beta_c_e)+0.25*(theta_ini-beta_c_e_ini)^2+0.25*(theta+beta_c_e)^2-alpha_c_e+(1-theta)*log2(1+Pr*(g2'*g2));

%     constraints(16)=alpha_c_e-alpha_c_1; % Linear scalar (real, 2 variables)
%     constraints(17)=alpha_c_e-alpha_c_2;

%     constraints(18)=(rho_c_1_ini^2)*(((h1*p_1)'*(h1*p_1))+((h1*p_2)'*(h1*p_2))+1)-2*real(p_c_ini'*h1'*h1*p_c)*rho_c_1_ini+((h1*p_c_ini)'*(h1*p_c_ini))*rho_c_1;
%     constraints(19)=(rho_c_2_ini^2)*(((h2*p_1)'*(h2*p_1))+((h2*p_2)'*(h2*p_2))+1)-2*real(p_c_ini'*h2'*h2*p_c)*rho_c_2_ini+((h2*p_c_ini)'*(h2*p_c_ini))*rho_c_2;
    constraints(18)=((h1*p_1)'*(h1*p_1))+((h1*p_2)'*(h1*p_2))+1-2*real(p_c_ini'*h1'*h1*p_c)/rho_c_1_ini+((h1*p_c_ini)'*(h1*p_c_ini))*rho_c_1/(rho_c_1_ini^2);
    constraints(19)=((h2*p_1)'*(h2*p_1))+((h2*p_2)'*(h2*p_2))+1-2*real(p_c_ini'*h2'*h2*p_c)/rho_c_2_ini+((h2*p_c_ini)'*(h2*p_c_ini))*rho_c_2/(rho_c_2_ini^2);
    constraints(20)=-(((g1*p_1_ini)'*(g1*p_1_ini)+(g1*p_2_ini)'*(g1*p_2_ini))*(rho_c_e-2*rho_c_e_ini)+rho_c_e+2*rho_c_e_ini*real(p_1_ini'*g1'*g1*p_1+p_2_ini'*g1'*g1*p_2))+(g1*p_c)'*(g1*p_c);

    constraints(21)=exp(beta_c_1*log(2))-(1+rho_c_1);
    constraints(22)=exp(beta_c_2*log(2))-(1+rho_c_2); 
    constraints(23)=1+rho_c_e-(1+log(2)*(beta_c_e-beta_c_e_ini))*exp(beta_c_e_ini*log(2)); 

    constraints(24)=-alpha_p_2;
    constraints(25)=-alpha_2_e;

    constraints(26)=-beta_p_2;
    constraints(27)=-beta_2_e;

    constraints(28)=-rho_p_2;
    constraints(29)=-rho_2_e;

    constraints(30)=-0.5*(theta_ini+beta_p_2_ini)*(theta+beta_p_2)+0.25*(theta_ini+beta_p_2_ini)^2+0.25*(theta-beta_p_2)^2+alpha_p_2;
    constraints(31)=-0.5*(theta_ini-beta_2_e_ini)*(theta-beta_2_e)+0.25*(theta_ini-beta_2_e_ini)^2+0.25*(theta+beta_2_e)^2-alpha_2_e;

%     constraints(32)=(rho_p_2_ini^2)*(((h2*p_1)'*(h2*p_1))+1)-2*real(p_2_ini'*h2'*h2*p_2)*rho_p_2_ini+((h2*p_2_ini)'*(h2*p_2_ini))*rho_p_2;  
    constraints(32)=((h2*p_1)'*(h2*p_1))+1-2*real(p_2_ini'*h2'*h2*p_2)/rho_p_2_ini+((h2*p_2_ini)'*(h2*p_2_ini))*rho_p_2/(rho_p_2_ini^2);  
    constraints(33)=-(((g1*p_c_ini)'*(g1*p_c_ini)+(g1*p_1_ini)'*(g1*p_1_ini))*(rho_2_e-2*rho_2_e_ini)+rho_2_e+2*rho_2_e_ini*real(p_c_ini'*g1'*g1*p_c+p_1_ini'*g1'*g1*p_1))+(g1*p_2)'*(g1*p_2);

    constraints(34)=exp(beta_p_2*log(2))-(1+rho_p_2);
    constraints(35)=1+rho_2_e-(1+log(2)*(beta_2_e-beta_2_e_ini))*exp(beta_2_e_ini*log(2));

    Constraints=[constraints(1:15)<=0,constraints(18:35)<=0];  

    P = optimize(Constraints,-Objective,ops);

    if abs(value(Objective)-t_past)<=tolerance
        loop=0;
    else
%         fprintf(' %4.0f : %12.3E \n',[count value(Objective)]);
        t_past=value(Objective);
%         disp_iter(1,count+1)=t_past;
        count=count+1;
        p_1_ini=value(p_1);
        p_2_ini=value(p_2);
        p_c_ini=value(p_c);
        theta_ini=value(theta);
        beta_p_1_ini=value(beta_p_1);
        beta_p_2_ini=value(beta_p_2);
        beta_1_e_ini=value(beta_1_e);
        beta_2_e_ini=value(beta_2_e);
        beta_c_1_ini=value(beta_c_1);
        beta_c_2_ini=value(beta_c_2);
        beta_c_e_ini=value(beta_c_e);

        rho_p_1_ini=value(rho_p_1);
        rho_p_2_ini=value(rho_p_2);
        rho_1_e_ini=value(rho_1_e);
        rho_2_e_ini=value(rho_2_e);
        rho_c_1_ini=value(rho_c_1);
        rho_c_2_ini=value(rho_c_2);
        rho_c_e_ini=value(rho_c_e);
    end
    
    if count>=num_iter
        break;
    end
   
end

%% show the outputs
% R_c=min(value(alpha_c_1),value(alpha_c_2));
% C_c_e=value(alpha_c_e);
% R_p_1=value(alpha_p_1);
% R_p_2=value(alpha_p_2);
% C_1_e=value(alpha_1_e);
% C_2_e=value(alpha_2_e);
ratio_comm=value(p_c'*p_c)/Pt;
theta_opt=value(theta);
SR_RS_3=value(Objective);
fprintf('RS : theta_3= %1.0f,rate_3= %1.3E,ratio_3=%1.3f \n',[theta_opt SR_RS_3 ratio_comm]);
end