function [SR_MU_LP_2]=sumRateMU_2(Pt,h1,h2,g1,NT,p_1_ini,p_2_ini,tolerance,num_iter)
%% introduced parameters
T_1=abs(h1*p_1_ini)^2+abs(h1*p_2_ini)^2+1;
T_2=abs(h2*p_1_ini)^2+abs(h2*p_2_ini)^2+1;
T_3=abs(g1*p_1_ini)^2+abs(g1*p_2_ini)^2+1;
T_4=abs(g1*p_1_ini)^2;
T_5=abs(g1*p_2_ini)^2;

I_1=T_1-abs(h1*p_1_ini)^2;
I_2=T_2-abs(h2*p_2_ini)^2;
I_3=T_3-abs(g1*p_1_ini)^2;
I_4=T_3-abs(g1*p_2_ini)^2;

% \rho_p=[rho_p_1_ini,rho_p_2_ini,rho_1_e_ini,rho_2_e_ini]
rho_p_1_ini=T_1/I_1-1;
rho_p_2_ini=T_2/I_2-1;
rho_1_e_ini=T_4/I_3;
rho_2_e_ini=T_5/I_4;

beta_1_e_ini=log2(1+rho_1_e_ini);
beta_2_e_ini=log2(1+rho_2_e_ini);


p_1=sdpvar(NT,1,'full','complex');
p_2=sdpvar(NT,1,'full','complex');

beta_p_1=sdpvar;
beta_p_2=sdpvar;
beta_1_e=sdpvar;
beta_2_e=sdpvar;   

rho_p_1=sdpvar;
rho_p_2=sdpvar;
rho_1_e=sdpvar;
rho_2_e=sdpvar; 

loop=1;
t_past=0;
count=0;

while (loop)     
    %constraints        
    constraints(1)=-beta_p_1;
    constraints(2)=-beta_p_2;
    constraints(3)=-beta_1_e;
    constraints(4)=-beta_2_e; 
    
    constraints(5)=-rho_p_1;
    constraints(6)=-rho_p_2;
    constraints(7)=-rho_1_e;
    constraints(8)=-rho_2_e;

    constraints(9)=p_1'*p_1+p_2'*p_2-Pt; % Quadratic scalar (real, 24 variables)
    
    % Quadratic scalar (real)
    constraints(10)=((h1*p_2)'*(h1*p_2))+1-2*real(p_1_ini'*h1'*h1*p_1)/rho_p_1_ini+((h1*p_1_ini)'*(h1*p_1_ini))*rho_p_1/(rho_p_1_ini^2);
    constraints(11)=((h2*p_1)'*(h2*p_1))+1-2*real(p_2_ini'*h2'*h2*p_2)/rho_p_2_ini+((h2*p_2_ini)'*(h2*p_2_ini))*rho_p_2/(rho_p_2_ini^2);  

    constraints(12)=-(((g1*p_2_ini)'*(g1*p_2_ini))*(rho_1_e-2*rho_1_e_ini)+rho_1_e+2*rho_1_e_ini*real(p_2_ini'*g1'*g1*p_2))+(g1*p_1)'*(g1*p_1);
    constraints(13)=-(((g1*p_1_ini)'*(g1*p_1_ini))*(rho_2_e-2*rho_2_e_ini)+rho_2_e+2*rho_2_e_ini*real(p_1_ini'*g1'*g1*p_1))+(g1*p_2)'*(g1*p_2);

    % Nonlinear scalar (real)
    constraints(14)=exp(beta_p_1*log(2))-(1+rho_p_1); 
    constraints(15)=exp(beta_p_2*log(2))-(1+rho_p_2);
    
    % Linear scalar (real)   
    constraints(16)=1+rho_1_e-(1+log(2)*(beta_1_e-beta_1_e_ini))*exp(beta_1_e_ini*log(2));
    constraints(17)=1+rho_2_e-(1+log(2)*(beta_2_e-beta_2_e_ini))*exp(beta_2_e_ini*log(2));  
    constraints(18)=beta_p_1-beta_1_e; 
    constraints(19)=-beta_p_2+beta_2_e; 
    Constraints=[constraints(1:19)<=0];
    
%     ops = sdpsettings('solver','fmincon','fmincon.Algorithm','sqp','fmincon.MaxIter',Inf,'fmincon.MaxFunEvals',Inf,'verbose',0);
    % ops = sdpsettings('solver','fmincon','fmincon.MaxIter',Inf,'fmincon.MaxFunEvals',Inf,'verbose',1);
    ops = sdpsettings('verbose',0);
    
    Objective=beta_p_1-beta_1_e+beta_p_2-beta_2_e; 
    
    P = optimize(Constraints,-Objective,ops);

    if abs(value(Objective)-t_past)<=tolerance
        loop=0;
    else
        t_past=value(Objective);
        count=count+1;
        
        p_1_ini=value(p_1);
        p_2_ini=value(p_2);
        beta_1_e_ini=value(beta_1_e);
        beta_2_e_ini=value(beta_2_e);

        rho_p_1_ini=value(rho_p_1);
        rho_p_2_ini=value(rho_p_2);
        rho_1_e_ini=value(rho_1_e);
        rho_2_e_ini=value(rho_2_e);
    end
    
    if count>=num_iter
        break;
    end
   
end

%% show the outputs
% R_p_1=value(beta_p_1);
% R_p_2=value(beta_p_2);
% C_1_e=value(beta_1_e);
% C_2_e=value(beta_2_e);
SR_MU_LP_2=value(Objective);
% pow_p_1=value(norm(p_1)^2);
% pow_p_2=value(norm(p_2)^2);
fprintf('MU_2 : rate=%1.3E \n',[SR_MU_LP_2]);
end