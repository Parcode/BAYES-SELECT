%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define model used, read experimental data. 
%% In this script start (MCSTART=0) equilibrate MCMC walker and 
%% generate fitting parameters in accordance to chi2
%% or with MCSTART=1, take previous run and generate fitting parameters in accordance to chi2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% READ experimental data:
READ_EXPERIMENT_v0c;
yerr=err_rel;
ydata=R1;
const.NN=30;
const.NNtypes=3;
gamma1H  = 267.513e6;
omega= B0.* 1.0e6*2*pi;
const.omega=omega;
const.kb=physconst('Boltzmann');
const.visc=[1306.9, 1138.2, 1002.0].*1e-6; %% Pa*s
const.beta_intra = 5.492036e10; 
const.Nfunc=0;
%%
%% Parametrization in v5:
%% tau_w1, tau_w2 => residence time in two water sites x(1),x(2)
%% f1*S_1, f2*S_2 => order parameters for two water sites    x(3), x(4)
%% alpha2 non-dispersive relaxation contributions     x(5)
%% R- micelle radius ( around 54Å)
%%
%% we need to define boundaries on fitted parameters:
xmin(1:3)   = 1.0e-5;    %% f1_S1
xmax(1:3)   = 1.0e-1;
xmin(4:6)   = 1.0e-5;    %% f2_S2
xmax(4:6)   = 1.0e-1;
xmin(7:9)   = 1.0e-9; %% tau_c1
xmax(7:9)   = 10e-6;
xmin(10:12) = 1e-10;  %% tau_c2
xmax(10:12) = 1e-7;
xmin(13:15) = 1e-3;   %% alpha2
xmax(13:15) = 0.80;
%%xmin(16)    = 0.1;   %% WD
%%xmax(16)    = 3.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MC loop                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nparm = length(xmax);
N_warmup = 1000;
N_MCMC   = 10000;
N_sub = 100;
MCSTART=0
NBATCH=0;
beta=[1.0];
%%beta=beta2;
const.beta=beta;
err=1e10.*ones(length(xmax),1);
errtot=1e10;
mcmcparm.err=err;
mcmcparm.errtot=errtot;
mcmcparm.yerr=yerr;
mcmcparm.N_warmup=N_warmup;
mcmcparm.N_MCMC=N_MCMC;
mcmcparm.N_sub=N_sub;
mcmcparm.MCSTART=MCSTART;
mcmcparm.xmin=xmin;
mcmcparm.xmax=xmax;
vars.mcmcparm=mcmcparm;
vars.const=const;
vars.x=(xmax-xmin);
vars.seed=715;
rng(vars.seed);
vars.ydata=ydata;
filen='MCMC_v0c_B/MCMC_v0c7_beta_';
%% GENERATE N_MCMC random numbers of form prior distribution (here uniform)
%%  [mcmcout]  = MCMCdriver_v1(vars);
X0 = xmin.*ones(N_MCMC,size(xmin,2)) + (xmax-xmin).*rand(N_MCMC,size(xmin,2));
for ii=1:N_MCMC,
  [errtot,err,ymodel] = Rfit_tot(X0(ii,:),ydata,const,yerr);
  err_MCMC(ii,1)=errtot;
end;
mcmcout.X0=X0;
mcmcout.err_MCMC=err_MCMC;
mcmcout.N_MCMC=N_MCMC;
mcmcout.const=const;
tmpname=strcat(filen,num2str(NBATCH),'.mat');
mcmcout.const.beta=0.0d0;
save(tmpname,'mcmcout');
return;
