%% plot exp data
READ_EXPERIMENT_v0c;
figure;
omega= B0.*2*pi*1e6;
semilogx(B0(:,1),R1(:,1),'kv','MarkerSize',12,'LineWidth',2);
%set(gca,'xscale','log');
hold on;
semilogx(B0(:,2),R1(:,2),'kx','MarkerSize',12,'LineWidth',2);
semilogx(B0(:,3),R1(:,3),'ko','MarkerSize',12,'LineWidth',2);
set(gca,'XScale','log');
set(gca,'fontsize', 24);
set(gca,'LineWidth',2);
xlabel('Resonance frequency (MHz)','FontSize',24);
ylabel('Relaxation rate (s^{-1})','FontSize',24)
%%lgd=legend('T=10 ^oC','T=15 ^oC','T=20 ^oC');
%%lgd.FontSize = 24;

%%
%%load  MCMC_v0c.mat;
%%load MCMC_v0c_beta_22.mat
%%Y0=X0;
%%err_tot=err_MCMC;
filetmp=cell(4,1);
tmpftmp{1}='MCMC_v0c_B/MCMC_v0c_beta_22.mat';
tmpftmp{2}='MCMC_v0c_B/MCMC_v0c2_beta_22.mat';
tmpftmp{3}='MCMC_v0c_B/MCMC_v0c3_beta_22.mat';
tmpftmp{4}='MCMC_v0c_B/MCMC_v0c4_beta_22.mat';
tmpftmp{5}='MCMC_v0c_B/MCMC_v0c5_beta_22.mat';
tmpftmp{6}='MCMC_v0c_B/MCMC_v0c6_beta_22.mat';
tmpftmp{7}='MCMC_v0c_B/MCMC_v0c7_beta_22.mat';
N_MCMC=0;
X0=[];
err_tot=[];
Yerrbest=1e10;
for ii=1:7,
  load(tmpftmp{ii});
  const=mcmcout.const;
  if(Yerrbest>mcmcout.errbest)
    Yerrbest=mcmcout.errbest
    Ybest=mcmcout.Xbest;
    ii
  end;
  N_MCMC=N_MCMC + mcmcout.N_MCMC;
  X0=[X0;mcmcout.X0];
  err_tot=[err_tot;mcmcout.err_MCMC];
end;
Y0=X0;
R1_R2model=@R1_R2model_v0c
const.Nfunc=0
yerr=err_rel;
ydata=R1;
const.NN=30;
const.NNtypes=3;
gamma1H  = 267.513e6;
const.omega=omega;
%%x=mean(Y0);
x=Ybest;
disp('f1_S1: ');
disp(x(1:3)*1e3);
disp('f2_S2: ');
disp(x(4:6)*1e3);
disp('tau_c1:');
disp(x(7:9)*1e9);
disp('tau_c2:');
disp(x(10:12)*1e9);
disp('alpha2:');
disp(x(13:15));

%%x=Xbest;
%%R1_R2model=const.R1_R2model;
 [errtot,err,ymodel] = Rfit_tot(x,ydata,const,yerr);
Nblock=10;
for ii=1:Nblock,
  nn=[(ii-1)*N_MCMC/Nblock+1,ii*N_MCMC/Nblock];
  chi2(ii)=2.0*mean(err_tot(nn(1):nn(2),1));
end;
fprintf('best chi2: %g  mean chi2: %g dev %g\n',errtot*2,mean(chi2),std(chi2)/sqrt(Nblock));

plot(B0(:,:),ymodel(:,:),'--r','linewidth',4);
str = ['A'];
%%tx = text(0.02,7.5,0,str);
tx = text(0.02,3,0,str);
tx.FontSize = 30;
tx.FontWeight = 'bold';
ylim([0 3.25])
set(gca,'fontsize', 24);
%%
for jj=1:const.NN,
    [ymodel(jj,1:const.NNtypes),RR1(jj,:),RR2(jj,:)] = R1_R2model(x,ydata,const,jj);
end;

figure;
semilogx(B0(:,3),ydata(:,3),'ko',B0(:,3),ymodel(:,3),'--r','LineWidth',4,'Markersize',12); 
hold on;
semilogx(B0(:,3),ones(30,1).*x(15),'-g','LineWidth',4);
semilogx(B0(:,3),RR1(:,3),'-b','LineWidth',4);
semilogx(B0(:,3),RR2(:,3),'-c','LineWidth',4);
str = ['B'];
tx = text(0.02,2.5,0,str);
tx.FontSize = 30;
tx.FontWeight = 'bold';
xlabel('Resonance frequency (MHz)','FontSize',20);
ylabel('Relaxation rate (s^{-1})','FontSize',20)
lgd=legend('EXP','Total Model','\alpha','R\beta{1}','R\beta{2}');
lgd.FontSize = 24;
set(gca,'LineWidth',2);
set(gca,'fontsize', 24);
ylim([0 3.0])
for nn=1:6,
  [f1,x1] = ecdf(Y0(:,nn).*1e3);
  iF=find(f1<0.005);
  xlow(nn)=x1(max(iF));
  iF=find(f1>1-0.005);
  xhigh(nn)=x1(min(iF));
  (xhigh(nn)-xlow(nn))/2.0;
  fprintf('%g  %g  %g  %g\n',x(nn)*1e3,(xhigh(nn)-xlow(nn))/2.0,xhigh(nn),xlow(nn));
end;

for nn=7:12,
  [f1,x1] = ecdf(Y0(:,nn).*1e9);
  iF=find(f1<0.005);
  xlow(nn)=x1(max(iF));
  iF=find(f1>1-0.005);
  xhigh(nn)=x1(min(iF));
  (xhigh(nn)-xlow(nn))/2.0;
  fprintf('%g  %g  %g  %g\n',x(nn)*1e9,(xhigh(nn)-xlow(nn))/2.0,xhigh(nn),xlow(nn));
end;

for nn=13:15,
[f1,x1] = ecdf(Y0(:,nn));
iF=find(f1<0.005);
xlow(nn)=x1(max(iF));
iF=find(f1>1-0.005);
xhigh(nn)=x1(min(iF));
(xhigh(nn)-xlow(nn))/2.0;
  fprintf('%g  %g  %g  %g\n',x(nn),(xhigh(nn)-xlow(nn))/2.0,xhigh(nn),xlow(nn));
end;

