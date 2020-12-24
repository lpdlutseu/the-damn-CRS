% color=[1 0 0;0 1 0;0 0 1;0.5 1 1;1 1 0.5;1 0.5 1; 0 0 0.5; 0.5 0 0;0 0.5 0;1 0.5 0.5];
% figure
% % plot(SNRdB,robustSR_SCA);
% 
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
MU_LP_SCA=mean(SR_MU_LP,2)';
RS_SCA=mean(SR_RS,2)';
plot(SNRdB,MU_LP_SCA,'-o',...
    SNRdB,RS_SCA,'-s','linewidth',2);
xlabel('SNR(dB)');
ylabel('Secrecy Sum Rate (bit/s/Hz)');
grid on;
legend('MU-LP','CRS','Location','southeast')
set(gca,'fontsize',12);
grid on;
print -deps epsFig