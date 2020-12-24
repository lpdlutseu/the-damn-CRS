% function [SR_RS_4,SR_MU_LP_3,SR_NOMA]= RS_SCA_rateRegion1(H,h3,g2,Pt,Pr,ind_relay,tolerance)
function [SR_RS_4]= RS_SCA_rateRegion1(H,h3,g2,Pt,Pr,ind_relay,tolerance,epsilon_1,epsilon_2)
%  function [SR_RS,SR_NOMA]= RS_SCA_rateRegion1(H,h3,g2,Pt,Pr,ind_relay,tolerance)
%% system parameters
[~,NT,N_user] = size(H) ;
h1=H(:,:,1); h2=H(:,:,2); g1=H(:,:,3); 

%% number of iterations
num_iter=600;

%% obtain the initial values of p_c_1_ini, p_c_2_ini, p_c_e_ini, theta_ini    
[p_1_ini,p_2_ini,p_c_ini,theta_ini]=seekforini(Pt,Pr,NT,N_user,h1,h2,h3,g1,g2,ind_relay);

%% RS-SCA
% SR_RS_1=sumRateRS_1(Pt,Pr,h1,h2,h3,g1,g2,NT,p_1_ini,p_2_ini,p_c_ini,theta_ini,ind_relay,tolerance,num_iter);
% SR_RS_2=sumRateRS_2(Pt,Pr,h1,h2,h3,g1,g2,NT,p_1_ini,p_2_ini,p_c_ini,theta_ini,ind_relay,tolerance,num_iter);
% SR_RS_3=sumRateRS_3(Pt,Pr,h1,h2,h3,g1,g2,NT,p_1_ini,p_2_ini,p_c_ini,theta_ini,ind_relay,tolerance,num_iter);
SR_RS_4=sumRateRS_4(Pt,Pr,h1,h2,h3,g1,g2,NT,p_1_ini,p_2_ini,p_c_ini,theta_ini,ind_relay,tolerance,num_iter,epsilon_1,epsilon_2);
% SR_RS=max([SR_RS_4]);
%% MU-LP-SCA
% [SR_MU_LP_1]=sumRateMU_1(Pt,h1,h2,g1,NT,p_1_ini,p_2_ini,tolerance,num_iter);
% [SR_MU_LP_2]=sumRateMU_2(Pt,h1,h2,g1,NT,p_1_ini,p_2_ini,tolerance,num_iter);
% [SR_MU_LP_3]=sumRateMU_3(Pt,h1,h2,g1,NT,p_1_ini,p_2_ini,tolerance,num_iter);
% SR_MU_LP=max([SR_MU_LP_1,SR_MU_LP_2,SR_MU_LP_3])

%% NOMA_SCA
% SR_NOMA_1=sumRateNOMA_1(Pt,Pr,h1,h2,h3,g1,g2,NT,p_1_ini,p_2_ini,p_c_ini,theta_ini,ind_relay,tolerance,num_iter);
% SR_NOMA_2=sumRateNOMA_2(Pt,Pr,h1,h2,h3,g1,g2,NT,p_1_ini,p_2_ini,p_c_ini,theta_ini,ind_relay,tolerance,num_iter);
% SR_NOMA=max([SR_NOMA_1,SR_NOMA_2]);

%% RS-SCA when theta=1,norm(g1*p_c)=0
% SR_RS=sumRateRS_new(Pt,Pr,h1,h2,h3,g1,g2,NT,p_1_ini,p_2_ini,p_c_ini,theta_ini,ind_relay,tolerance,num_iter);
%% RS-NRS theta=1
% SR_NRS=sumRateNRS(Pt,h1,h2,g1,NT,p_1_ini,p_2_ini,p_c_ini,tolerance,num_iter)
end
