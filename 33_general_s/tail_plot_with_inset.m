close all;
logt=log10(output.time_arr);
psi_extr=zeros(size(output.time_arr));
[r,c]=find_value_arr(dg_globals.r_sch,100)
% r=dg_globals.Np; c=dg_globals.K;
for i=1:length(output.time_arr)
   psi_t=flatten_array(output.psi_arr(i,:,:),dg_globals);
  
   psi=psi_t(r(1),c(1));


   psi_extr(i)=psi;
end


% subplot(2,1,1)
plot(logt,log10(psi_extr),"m.",LineWidth=2); 
xlim([1.4,3.3])
title("s=0")
% ylim([-15, 2 ])

% subplot(2,1,2)
% 
% log_psi=log10(psi_extr);
% plot(logt,log10(psi_extr),"m-",LineWidth=1.5); 
% xlim([1.4,3.3])
% ylim([-15, 2 ])
% 
% hold on

% % output_s1=load('s1.mat');
% psi_extr_s1=zeros(size(output_s1.time_arr));
% logt_s1=log10(output_s1.time_arr);
% [r,c]=find_value_arr(dg_globals.r_sch,100)
% % r=dg_globals.Np; c=dg_globals.K;
% for i=1:length(output_s1.time_arr)
%    psi_t_s1=flatten_array(output_s1.psi_arr(i,:,:),dg_globals);
% 
%    psi_s1=psi_t_s1(r(1),c(1));
% 
% 
%    psi_extr_s1(i)=psi_s1;
% end
% 
% log_psi_s1=log10(psi_extr_s1);
% plot(logt_s1,log10(psi_extr_s1),"b-.",LineWidth=1.5); 
% 
% 
% 
% 
% % output_s2=load('s2.mat');
% psi_extr_s2=zeros(size(output_s2.time_arr));
% logt_s2=log10(output_s2.time_arr);
% [r,c]=find_value_arr(dg_globals.r_sch,100)
% % r=dg_globals.Np; c=dg_globals.K;
% for i=1:length(output_s2.time_arr)
%    psi_t_s2=flatten_array(output_s2.psi_arr(i,:,:),dg_globals);
% 
%    psi_s2=psi_t_s2(r(1),c(1));
% 
% 
%    psi_extr_s2(i)=psi_s2;
% end
% 
% log_psi_s2=log10(psi_extr_s2);
% plot(logt_s2,log10(psi_extr_s2),"k-.",LineWidth=1.5); 
% hold off
% 
% 
% legend("s=0","s=-1","s=-2")
% legend('Location','best')


% cut1=find(abs(logt-3)<1e-3);
% cut1=cut1(end);
% cut2=find(abs(logt-3.9)<1e-3);
% cut2=cut2(end);
% % subplot(2,1,1)
% 
% plot(logt(4:end),log10(abs((psi_extr0(4:end)))),"r-",'LineWidth',1.5);
dLog=ComputeTailDecayRate(output.time_arr,psi_extr);
% % plot(logt(cut1:cut2),dLog0(cut1:cut2),'r','LineWidth',1.5)
% data_fit_x=output.tarr(cut1:end-1)'; data_fit_y=dLog0(cut1:end)';
% P0=fit(data_fit_x(:),data_fit_y(:),"c+b*x^(-m)");
% expected_rate0=P0(1e6);
% hold on 
% % plot(logt(cut1:end-1),log10(abs((psi_extr0(cut1:end-1)))),"r-*",'LineWidth',1.5)
% % hold on
% 
% 
% plot(logt(4:end),log10(abs((psi_extr1(4:end)))),"k--",'LineWidth',1.5);
% hold on
% dLog1=ComputeTailDecayRate(output.tarr,psi_extr1);
% data_fit_x=output.tarr(cut1:end-1)'; data_fit_y=dLog1(cut1:end)';
% P1=fit(data_fit_x(:),data_fit_y(:),"c+b*x^(-m)");
% expected_rate1=P1(1e6);
% 
% plot(logt(4:end) ,log10(abs((psi_extr2(4:end)))),"b-.",'LineWidth',1.5);
% hold on
% dLog2=ComputeTailDecayRate(output.tarr,psi_extr2);
% data_fit_x=output.tarr(cut1:end-1)'; data_fit_y=dLog2(cut1:end)';
% P2=fit(data_fit_x(:),data_fit_y(:),"c+b*x^(-m)");
% expected_rate2=P2(1e6);
% % xline(logt(cut1))
% % xline(logt(cut2))
% % plot(logt(cut1:end-1),log10(abs((psi_extr2(cut1:end-1)))),"b--",'LineWidth',1.5)
% % hold on
% 
% ylabel("$\log_{10} |\Psi (\tau,\rho)|$",'interpreter','latex','fontsize',14);
% xlabel("$\log_{10}(\tau)$",'interpreter','latex','fontsize',14)
% % t=sprintf("N=%d, subd=%d, FinalTime=%d, \nsignal extracted at $$r_{*}=%d $$ \nwith moving initial data",dg_globals.N,dg_globals.subd,FinalTime,dg_globals.r_adj(r(1),c(1)));
% % title(sprintf("Schwarzschild Tail for \nNon-Zero Momentum Initial Conditon \n at $$r_{*}=%1.2f M $$",dg_globals.r_adj(r(1),c(1))),'interpreter','latex')
% % legend_str=sprintf("$\Psi_{0} (%1.3f)$","$\Psi_{1}$","$\Psi_{2}$",'interpreter','latex')
% legend("$\Psi_{00}$" + sprintf("(-3,%1.6f)",expected_rate0),"$\Psi_{10}$"+ sprintf("(-5.%1.6f)",expected_rate1),"$\Psi_{20}$"+ sprintf("(-7.%1.6f)",expected_rate2),'interpreter','latex')
% legend('Location','best')
% grid()
% ylim([-20 1])
% set(gca,'fontsize',12)
% 
% 
% 
% % % % create smaller axes in top right, and plot on it
axes('Position',[.26 .2 .47 .2])
box on
grid()
cut1_inplot=find(abs(logt-2.4)<1e-3);
cut1_inplot=cut1_inplot(1);
cut2_inplot=find(abs(logt-3)<1e-3);
cut2_inplot=cut2_inplot(end);
% plot(logt(cut1_inplot:cut2_inplot),dLog(cut1_inplot:cut2_inplot),"r-",'LineWidth',1.5);
plot(logt(2:end),dLog,"r-",'LineWidth',1.5);
grid()
xlim([2.4,2.9])

% yline(-3,"--")
% hold on
% plot(logt(cut1_inplot:cut2_inplot),dLog1(cut1_inplot:cut2_inplot),"k--",'LineWidth',1.5);
% yline(-5,"--")
% plot(logt(cut1_inplot:cut2_inplot),dLog2(cut1_inplot:cut2_inplot),"b-.",'LineWidth',1.5);
% yline(-7,"--")
% ylim([-9 -1])
% xlim([2.5 3.9])
% ylabel("LPI",'interpreter','latex')
% hold off

% ps=get(gcf, 'Position');
% save_fig_and_crop(ps)