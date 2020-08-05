load MCMC_v0c_B/beta_v0c_beta_1.dat
load MCMC_v0c_B/beta_v0c_beta_2.dat
load MCMC_v0c_B/beta_v0c_beta_3.dat
load MCMC_v0c_B/beta_v0c_beta_4.dat
load MCMC_v0c_B/beta_v0c_beta_5.dat
load MCMC_v0c_B/beta_v0c_beta_6.dat
load MCMC_v0c_B/beta_v0c_beta_7.dat
Nfiles=7;
Ndata=28;
outputs=zeros(Ndata,Nfiles);
outputs(:,1)=beta_v0c_beta_1(:,2);
outputs(:,2)=beta_v0c_beta_3(:,2);
outputs(:,3)=beta_v0c_beta_4(:,2);
outputs(:,4)=beta_v0c_beta_5(:,2);
outputs(:,5)=beta_v0c_beta_7(:,2);
outputs(:,7)=beta_v0c_beta_6(:,2);
outputs(:,6)=beta_v0c_beta_2(:,2);
integ=mean(outputs,2);
stdint=(std(outputs')')./sqrt(Nfiles);
stdcumtrapz0=std(cumtrapz(beta_v0c_beta_1(:,1),outputs)')';
cumtrapz0=cumtrapz(beta_v0c_beta_1(:,1),integ);
Nfiles10=7;
load MCMC_v10c_B/beta_v10c_beta_1.dat
load MCMC_v10c_B/beta_v10c_beta_2.dat
load MCMC_v10c_B/beta_v10c_beta_3.dat
load MCMC_v10c_B/beta_v10c_beta_4.dat
load MCMC_v10c_B/beta_v10c_beta_5.dat
load MCMC_v10c_B/beta_v10c_beta_6.dat
load MCMC_v10c_B/beta_v10c_beta_7.dat
outputs10=zeros(Ndata,Nfiles10);
outputs10(:,1)=beta_v10c_beta_1(:,2);
outputs10(:,2)=beta_v10c_beta_2(:,2);
outputs10(:,3)=beta_v10c_beta_3(:,2);
outputs10(:,4)=beta_v10c_beta_4(:,2);
outputs10(:,5)=beta_v10c_beta_5(:,2);
outputs10(:,6)=beta_v10c_beta_6(:,2);
outputs10(:,7)=beta_v10c_beta_7(:,2);
integ10=mean(outputs10,2);
stdint10=(std(outputs10')')./sqrt(Nfiles10);
stdcumtrapz10=std(cumtrapz(beta_v10c_beta_1(:,1),outputs10)')';
cumtrapz10=cumtrapz(beta_v10c_beta_1(:,1),integ10);
figure;
errorbar(beta_v0c_beta_1(:,1),integ,stdint,'--or','linewidth',3);
hold on;
errorbar(beta_v10c_beta_1(:,1),integ10,stdint10,'--xk','linewidth',3);
set(gca,'yScale','log');
set(gca,'xScale','log');
xlabel('\beta')
legend('$\langle\frac{\chi^2_1({\bf x})}{2}\rangle_{\beta}$','$\langle\frac{\chi^2_3({\bf x})}{2}\rangle_{\beta}$','fontsize', ...
24,'Interpreter','latex','Location','southwest');
set(gca,'fontsize', 24);
str = ['A'];
ylim([10 1e12])
%%tx = text(1e-9,1e7,0,str);
%%tx.FontSize = 30;
%%tx.FontWeight = 'bold';
axes('Position',[0.58 0.6 0.3 0.3]);
box on;
errorbar(beta_v0c_beta_1(23:28,1),integ(23:28),stdint(23:28),'--or','linewidth',3);
hold on;
errorbar(beta_v10c_beta_1(23:28,1),integ10(23:28),stdint10(23:28),'--xk','linewidth',3);
ylim([10 70]);
xlim([0.2 1.1]);
set(gca,'fontsize', 24);

figure;
errorbar(beta_v0c_beta_1(:,1),cumtrapz0(:),stdcumtrapz0(:),'--or','linewidth',3);
hold on;
errorbar(beta_v10c_beta_1(:,1),cumtrapz10,stdcumtrapz10(:),'--xk','linewidth',3);
set(gca,'xScale','log');
xlabel('\beta')
legend('$\ln[Z_1(\beta)]$','$\ln[Z_3(\beta)]$','fontsize', 24,'Interpreter','latex','Location','northwest');
set(gca,'fontsize', 24);
str = ['B'];
%%tx = text(0.85,85,0,str);
%%tx.FontSize = 30;
%%tx.FontWeight = 'bold';

axes('Position',[0.25 0.45 0.3 0.3]);
box on;
errorbar(beta_v0c_beta_1(:,1),cumtrapz0(:),stdcumtrapz0(:),'--or','linewidth',3);
hold on;
errorbar(beta_v10c_beta_1(:,1),cumtrapz10,stdcumtrapz10(:),'--xk','linewidth',3);
ylim([1 350]);
xlim([0.01 1.1]);
set(gca,'fontsize', 24);

lz0=zeros(Nfiles,1);
lz10=zeros(3,1);

for ii=1:Nfiles,
  lz0(ii)=trapz(beta_v0c_beta_1(:,1),outputs(:,ii))
end;
for ii=1:Nfiles10,
  lz10(ii)=trapz(beta_v10c_beta_1(:,1),outputs10(:,ii))
end;

fprintf('lnZ_0c = %g +/- %g\n',mean(lz0),std(lz0)/sqrt(Nfiles));
fprintf('lnZ_10c = %g +/- %g\n',mean(lz10),std(lz10)/sqrt(Nfiles10));
