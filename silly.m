function SR_RS=sumRateRS(Pt,Pr,h1,h2,h3,g1,g2,NT,p_1_ini,p_2_ini,p_c_ini,theta_ini,ind_relay,tolerance,num_iter)
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
%% sdpvart defination
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

%% circulations and iterations
loop=ones(1,4);
t_past=zeros(1,4);

theta_ini_array=ones(1,4)*theta_ini;
p_c_ini_array=p_c_ini*ones(1,4);
p_1_ini_array=p_1_ini*ones(1,4);
p_2_ini_array=p_2_ini*ones(1,4);

% \rho_c=[rho_c_1_ini,rho_c_2_ini,rho_c_e_ini]
rho_c_1_ini_array=ones(1,4)*rho_c_1_ini;
rho_c_2_ini_array=ones(1,4)*rho_c_2_ini;
rho_c_e_ini_array=ones(1,4)*rho_c_e_ini;

% \rho_p=[rho_p_1_ini,rho_p_2_ini,rho_1_e_ini,rho_2_e_ini]
rho_p_1_ini_array=ones(1,4)*rho_p_1_ini;
rho_p_2_ini_array=ones(1,4)*rho_p_2_ini;
rho_1_e_ini_array=ones(1,4)*rho_1_e_ini;
rho_2_e_ini_array=ones(1,4)*rho_2_e_ini;

beta_c_1_ini_array=ones(1,4)*beta_c_1_ini;
beta_c_2_ini_array=ones(1,4)*beta_c_2_ini;
beta_c_e_ini_array=ones(1,4)*beta_c_e_ini;

beta_p_1_ini_array=ones(1,4)*beta_p_1_ini;
beta_p_2_ini_array=ones(1,4)*beta_p_2_ini;
beta_1_e_ini_array=ones(1,4)*beta_1_e_ini;
beta_2_e_ini_array=ones(1,4)*beta_2_e_ini;

Objective4=min(alpha_c_1,alpha_c_2)-alpha_c_e+alpha_p_1-alpha_1_e+alpha_p_2-alpha_2_e; 
Objective2=min(alpha_c_1,alpha_c_2)-alpha_c_e+alpha_p_1-alpha_1_e; 
Objective3=min(alpha_c_1,alpha_c_2)-alpha_c_e+alpha_p_2-alpha_2_e; 
Objective1=min(alpha_c_1,alpha_c_2)-alpha_c_e;
value_objective=[value(Objective1),value(Objective2),value(Objective3),value(Objective4)];

count=zeros(1,4);

for i=1:length(value_objective)
    while (loop(i))
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
            constraints(13)=-0.5*(theta_ini_array(i)+beta_c_1_ini_array(i))*(theta+beta_c_1)+0.25*(theta_ini_array(i)+beta_c_1_ini_array(i))^2+0.25*(theta-beta_c_1)^2+alpha_c_1;
            constraints(14)=-0.5*(theta_ini_array(i)+beta_c_2_ini_array(i))*(theta+beta_c_2)+0.25*(theta_ini_array(i)+beta_c_2_ini_array(i))^2+0.25*(theta-beta_c_2)^2-(1-theta)*log2(1+Pr*(h3'*h3))+alpha_c_2;  
         else
            constraints(13)=-0.5*(theta_ini_array(i)+beta_c_1_ini_array(i))*(theta+beta_c_1)+0.25*(theta_ini_array(i)+beta_c_1_ini_array(i))^2+0.25*(theta-beta_c_1)^2-(1-theta)*log2(1+Pr*(h3'*h3))+alpha_c_1;
            constraints(14)=-0.5*(theta_ini_array(i)+beta_c_2_ini_array(i))*(theta+beta_c_2)+0.25*(theta_ini_array(i)+beta_c_2_ini_array(i))^2+0.25*(theta-beta_c_2)^2+alpha_c_2;
        end

        constraints(15)=-0.5*(theta_ini_array(i)-beta_c_e_ini_array(i))*(theta-beta_c_e)+0.25*(theta_ini_array(i)-beta_c_e_ini_array(i))^2+0.25*(theta+beta_c_e)^2-alpha_c_e+(1-theta)*log2(1+Pr*(g2'*g2));

        constraints(16)=alpha_c_e-alpha_c_1; % Linear scalar (real, 2 variables)
        constraints(17)=alpha_c_e-alpha_c_2;

        constraints(18)=(rho_c_1_ini_array(i)^2)*(((h1*p_1)'*(h1*p_1))+((h1*p_2)'*(h1*p_2))+1)-2*real(p_c_ini_array(:,i)'*h1'*h1*p_c)*rho_c_1_ini_array(i)+((h1*p_c_ini_array(:,i))'*(h1*p_c_ini_array(:,i)))*rho_c_1;
        constraints(19)=(rho_c_2_ini_array(i)^2)*(((h2*p_1)'*(h2*p_1))+((h2*p_2)'*(h2*p_2))+1)-2*real(p_c_ini_array(:,i)'*h2'*h2*p_c)*rho_c_2_ini_array(i)+((h2*p_c_ini_array(:,i))'*(h2*p_c_ini_array(:,i)))*rho_c_2;
        constraints(20)=-(((g1*p_1_ini_array(:,i))'*(g1*p_1_ini_array(:,i))+(g1*p_2_ini_array(:,i))'*(g1*p_2_ini_array(:,i)))*(rho_c_e-2*rho_c_e_ini_array(i))+rho_c_e+2*rho_c_e_ini_array(i)*real(p_1_ini_array(:,i)'*g1'*g1*p_1+p_2_ini_array(:,i)'*g1'*g1*p_2))+(g1*p_c)'*(g1*p_c);

        constraints(21)=exp(beta_c_1*log(2))-(1+rho_c_1);
        constraints(22)=exp(beta_c_2*log(2))-(1+rho_c_2); 
        constraints(23)=1+rho_c_e-(1+log(2)*(beta_c_e-beta_c_e_ini_array(i)))*exp(beta_c_e_ini_array(i)*log(2)); 

        constraints(24)=-alpha_p_1;
        constraints(25)=-alpha_1_e;

        constraints(26)=-beta_p_1;
        constraints(27)=-beta_1_e;

        constraints(28)=-rho_p_1;
        constraints(29)=-rho_1_e;

        constraints(30)=-0.5*(theta_ini_array(i)+beta_p_1_ini_array(i))*(theta+beta_p_1)+0.25*(theta_ini_array(i)+beta_p_1_ini_array(i))^2+0.25*(theta-beta_p_1)^2+alpha_p_1; % Quadratic scalar (real, 3 variables)
        constraints(31)=-0.5*(theta_ini_array(i)-beta_1_e_ini_array(i))*(theta-beta_1_e)+0.25*(theta_ini_array(i)-beta_1_e_ini_array(i))^2+0.25*(theta+beta_1_e)^2-alpha_1_e;

        % Quadratic scalar (real) and Signomial scalar(real)
        constraints(32)=(rho_p_1_ini_array(i)^2)*(((h1*p_2)'*(h1*p_2))+1)-2*real(p_1_ini_array(:,i)'*h1'*h1*p_1)*rho_p_1_ini_array(i)+((h1*p_1_ini_array(:,i))'*(h1*p_1_ini_array(:,i)))*rho_p_1;
        constraints(33)=-(((g1*p_c_ini_array(:,i))'*(g1*p_c_ini_array(:,i))+(g1*p_2_ini_array(:,i))'*(g1*p_2_ini_array(:,i)))*(rho_1_e-2*rho_1_e_ini_array(i))+rho_1_e+2*rho_1_e_ini_array(i)*real(p_c_ini_array(:,i)'*g1'*g1*p_c+p_2_ini_array(:,i)'*g1'*g1*p_2))+(g1*p_1)'*(g1*p_1);

        constraints(34)=exp(beta_p_1*log(2))-(1+rho_p_1);    % Nonlinear scalar (real)   
        constraints(35)=1+rho_1_e-(1+log(2)*(beta_1_e-beta_1_e_ini_array(i)))*exp(beta_1_e_ini_array(i)*log(2));    % Linear scalar (real) 

        constraints(36)=-alpha_p_2;
        constraints(37)=-alpha_2_e;

        constraints(38)=-beta_p_2;
        constraints(39)=-beta_2_e;

        constraints(40)=-rho_p_2;
        constraints(41)=-rho_2_e;

        constraints(42)=-0.5*(theta_ini_array(i)+beta_p_2_ini_array(i))*(theta+beta_p_2)+0.25*(theta_ini_array(i)+beta_p_2_ini_array(i))^2+0.25*(theta-beta_p_2)^2+alpha_p_2;
        constraints(43)=-0.5*(theta_ini_array(i)-beta_2_e_ini_array(i))*(theta-beta_2_e)+0.25*(theta_ini_array(i)-beta_2_e_ini_array(i))^2+0.25*(theta+beta_2_e)^2-alpha_2_e;

        constraints(44)=(rho_p_2_ini_array(i)^2)*(((h2*p_1)'*(h2*p_1))+1)-2*real(p_2_ini_array(:,i)'*h2'*h2*p_2)*rho_p_2_ini_array(i)+((h2*p_2_ini_array(:,i))'*(h2*p_2_ini_array(:,i)))*rho_p_2;  
        constraints(45)=-(((g1*p_c_ini_array(:,i))'*(g1*p_c_ini_array(:,i))+(g1*p_1_ini_array(:,i))'*(g1*p_1_ini_array(:,i)))*(rho_2_e-2*rho_2_e_ini_array(i))+rho_2_e+2*rho_2_e_ini_array(i)*real(p_c_ini_array(:,i)'*g1'*g1*p_c+p_1_ini_array(:,i)'*g1'*g1*p_1))+(g1*p_2)'*(g1*p_2);

        constraints(46)=exp(beta_p_2*log(2))-(1+rho_p_2);
        constraints(47)=1+rho_2_e-(1+log(2)*(beta_2_e-beta_2_e_ini_array(i)))*exp(beta_2_e_ini_array(i)*log(2));
        
        if i==1
            Constraints1=[constraints(1:23)<=0]; %alpha_c_1,alpha_c_2,alpha_c_e
            P1 = optimize(Constraints1,-Objective1,ops);
        end
        
        if i==2
            Constraints2=[constraints(1:35)<=0]; %alpha_c_1,alpha_c_2,alpha_c_e,alpha_p_1,alpha_1_e
            P2 = optimize(Constraints2,-Objective2,ops);
        end
        
        if i==3
            Constraints3=[constraints(1:23,36:47)<=0]; %alpha_c_1,alpha_c_2,alpha_c_e,alpha_p_2,alpha_2_e
            P3 = optimize(Constraints3,-Objective3,ops);
        end
        
        if i==4
            Constraints4=[constraints(1:47)<=0]; %alpha_c_1,alpha_c_2,alpha_c_e,alpha_p_1,alpha_1_e,alpha_p_2,alpha_2_e
            P4 = optimize(Constraints4,-Objective4,ops);
        end

        if abs(value_objective(i)-t_past(i))<=tolerance
            loop(i)=0;
        else
            t_past(i)=value_objective(i);
            count(i)=count(i)+1;

            p_1_ini_array(:,i)=value(p_1);
            p_2_ini_array(:,i)=value(p_2);
            p_c_ini_array(:,i)=value(p_c);
            theta_ini_array(i)=value(theta);
            beta_p_1_ini_array(i)=value(beta_p_1);
            beta_p_2_ini_array(i)=value(beta_p_2);
            beta_1_e_ini_array(i)=value(beta_1_e);
            beta_2_e_ini_array(i)=value(beta_2_e);
            beta_c_1_ini_array(i)=value(beta_c_1);
            beta_c_2_ini_array(i)=value(beta_c_2);
            beta_c_e_ini_array(i)=value(beta_c_e);

            rho_p_1_ini_array(i)=value(rho_p_1);
            rho_p_2_ini_array(i)=value(rho_p_2);
            rho_1_e_ini_array(i)=value(rho_1_e);
            rho_2_e_ini_array(i)=value(rho_2_e);
            rho_c_1_ini_array(i)=value(rho_c_1);
            rho_c_2_ini_array(i)=value(rho_c_2);
            rho_c_e_ini_array(i)=value(rho_c_e);
        end

        if count(i)>=num_iter
            break;
        end
    end
   
end

%% show the outputs
% R_c=min(value(alpha_c_1),value(alpha_c_2));
% C_c_e=value(alpha_c_e);
% R_p_1=value(alpha_p_1);
% R_p_2=value(alpha_p_2);
% C_1_e=value(alpha_1_e);
% C_2_e=value(alpha_2_e);
% ratio_comm=value(p_c'*p_c)/Pt;
% theta_opt=value(theta);
SR_RS=max(value_objective);
fprintf('RS : theta= %1.0f,rate= %1.3E,ratio=%1.3f \n',[theta_opt SR_RS ratio_comm]);
end