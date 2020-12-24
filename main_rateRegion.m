 clc;clear all; clf
tic
%% parameter settings
NT=4; N_user=2; A_U=1;

tolerance = 10^-3;

% delta_g^H*delta_g<=epsilon_1^2
epsilon_1 = 0.01; epsilon_2 = 0.01;

bias_1=1; bias_2=1; bias_3=1; bias_e_1=0.1; bias_e_2=0.1;

SNRdB=[5:5:30];

N_channel=1;

SR_RS=zeros(length(SNRdB),N_channel);
% SR_MU_LP=zeros(length(SNRdB),N_channel);
% SR_NOMA=zeros(length(SNRdB),N_channel);

%% circulation
for i_channel=1:N_channel
   randn('seed',(i_channel)*3*N_user) 
   %% channel realization
   %Phase I
   H_BC(:,:,1)=sqrt(bias_1)/sqrt(2)*(randn(1,NT)+1i*randn(1,NT)); %h1
   H_BC(:,:,2)=sqrt(bias_2)/sqrt(2)*(randn(1,NT)+1i*randn(1,NT)); %h2
   H_BC(:,:,3)=sqrt(bias_e_1)/sqrt(2)*(randn(1,NT)+1i*randn(1,NT)); %g1,g_1^{^},the estimation of g1

   %Phase II
   h3=norm(sqrt(bias_3)/sqrt(2)*(randn(1,NT)+1i*randn(1,NT))); %h3
   g2=norm(sqrt(bias_e_2)/sqrt(2)*(randn(1,NT)+1i*randn(1,NT))); %g2,g_2^{^},the estimation of g2

   %% Relay selection
    if norm(H_BC(:,:,1)) >= norm(H_BC(:,:,2))
        ind_relay=1;
        fprintf('ind_relay=%1.0f \n',ind_relay);
    else
        ind_relay=2;
        fprintf('ind_relay=%1.0f \n',ind_relay);
    end
       
   for i_snr=1:length(SNRdB)
      %% transmit power  
       Pt=10^(SNRdB(i_snr)/10);
       Pr=Pt;
       fprintf('i_channel=%1.0f,i_snr=%1.0f \n',[i_channel,i_snr]);

        %% use yalmip
        %[SR_RS(i_snr,i_channel),SR_MU_LP(i_snr,i_channel),SR_NOMA(i_snr,i_channel)]= RS_SCA_rateRegion1(H_BC,h3,g2,Pt,Pr,ind_relay,tolerance);
        [SR_RS(i_snr,i_channel)]= RS_SCA_rateRegion1(H_BC,h3,g2,Pt,Pr,ind_relay,tolerance,epsilon_1,epsilon_2);
   end 
    
end 
toc
% save('data_1_1.mat','SR_RS','SR_MU_LP','SR_NOMA')
%% plot the average sum-rate
RS_SCA=mean(SR_RS,2)';
% MU_LP_SCA=mean(SR_MU_LP,2)';
% NOMA_SCA=mean(SR_NOMA,2)';
plot(SNRdB,RS_SCA,'-o');
% plot(SNRdB,RS_SCA,'-o',...
%     SNRdB,MU_LP_SCA,'-s',...
%     SNRdB,NOMA_SCA,'-*','linewidth',2);
xlabel('SNR(dB)');
ylabel('Secrecy Sum Rate (bit/s/Hz)');
grid on;
% legend('CRS','MU-LP','C-NOMA','Location','southeast')
set(gca,'fontsize',12);
grid on;
print -deps epsFig

% plot(SNRdB,SR_MU_LP,'-s','linewidth',2);
% MU_LP_SCA=mean(SR_MU_LP,2)';
% plot(SNRdB,MU_LP_SCA,'-s','linewidth',2);
%% plot each channel
% color=[1 0 0;0 1 0;0 0 1;0.5 1 1;1 1 0.5;1 0.5 1; 0 0 0.5; 0.5 0 0;0 0.5 0;1 0.5 0.5];
% figure
% plot(SNRdB,robustSR_SCA);

% for i_channel=1:N_channel    
%     plot(SNRdB,SR_MU_LP(:,i_channel)','color',color(i_channel,:));    
%     leg_str{i_channel}=['channel',num2str(i_channel),];    
%     hold on
% end
% 
% for i_channel=1:N_channel    
%     plot(SNRdB,SR_RS(:,i_channel)','color',color(i_channel,:));    
%     leg_str{i_channel}=['channel',num2str(i_channel),];    
%     hold on
% end
% legend(leg_str)
% 
% xlabel('SNR(dB)');
% ylabel('Secrecy Sum Rate (bit/s/Hz)');
% grid on;