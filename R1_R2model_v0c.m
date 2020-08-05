%%
%% R1_R2model_v1, use previous NMRD fortran formula (DD-inter & DD-intra)
%%
function [ymodel,RR1,RR2] = R1_R2model_v0c(x,ydata,const,jj)
%%
  w0=const.omega(jj,:); %% w0=omega*2.0d0*pi*1.0e6
%% Parametrization in v1:
%% tau_w1, tau_w2 => residence time in two water sites x(1),x(2)
%% f1*S_1, f2*S_2 => order parameters for two water sites    x(3), x(4)
%% alpha2 non-dispersive relaxation contributions     x(5)
%% R- micelle radius ( around 54Å)
%%
%% we need to define boundaries on fitted parameters:
%% xmin(1:3)   = 0.0;    %% f1_S1
%%(1:3)   = 1.0e-1;
%%(4:6)   = 0.0;    %% f2_S2
%%(4:6)   = 1.0e-1;
%%(7:9)   = 10e-12; %% tau_c1
%%(7:9)   = 10e-6;
%%(10:12) = 1e-12;  %% tau_c2
%%(10:12) = 10e-9;
%%(13:15) = 1e-3;   %% alpha2
%%(13:15) = 0.65;
%%(16)    = 0.05;   %% WD
%%(16)    = 0.5;
%%  R3=x(6)*1e-30;
%%  T=[10,15,20] + 273.15;
%%  tauR=((4*pi)/3.0)*const.visc.*R3./(const.kb.*T);
%%  tauC2= x(8).*tauR./(tauR+x(8));
  %%tauC =  x(9).*tauR./(tauR+x(9));%%[x(9),x(10),x(11)];
  f1_S1=x(1:3);
  f2_S2=x(4:6);
  tauC=x(7:9);
  tauC2=x(10:12);
  alpha2=x(13:15);
  WD=0.4;
  beta_intra=const.beta_intra;
%%%%
  omw2=w0.*w0.*tauC.*tauC;
  omw22=w0.*w0.*tauC2.*tauC2;
%% "const.NNtypes" number of temperatures per resonance frequency:
%%  RR1intra(1,1:const.NNtypes)=(0.2./(1+omw2(1,1:const.NNtypes))+0.8./(1+4.*omw2(1,1:const.NNtypes)));
%%  RR1inter(1,1:const.NNtypes)=(0.1+0.3./(1+omw2(1,1:const.NNtypes))+0.6./(1+4.*omw2(1,1:const.NNtypes)));
%%  RR1ampl(1,1:const.NNtypes)=(f1_S1(1:const.NNtypes)*beta_intra.*tauC(1,1:const.NNtypes)-alpha2(1,1:const.NNtypes))./(1.0+WD);
%%  RR1(1,1:const.NNtypes)= RR1ampl(1,:).*( RR1intra(1,:) + WD*RR1inter(1,:) );
%%
%%  RR2intra=(0.2./(1+omw22(1,1:const.NNtypes))+0.8./(1+4.*omw22(1,1:const.NNtypes)));
%%  RR2inter=(0.1+0.3./(1+omw22(1,1:const.NNtypes))+0.6./(1+4.*omw22(1,1:const.NNtypes)));
%%  RR2ampl(1,1:const.NNtypes)=(f2_S2(1:const.NNtypes).*beta_intra.*tauC2(1,1:const.NNtypes)-alpha2(1,1:const.NNtypes))./(1.0+WD);
%%  RR2(1,1:const.NNtypes)= (RR2ampl(1,1:const.NNtypes)).*( RR2intra(1,1:const.NNtypes) + WD.*RR2inter(1,1:const.NNtypes ));
  for ii=1:const.NNtypes,
    RR1intra=(0.2/(1+omw2(1,ii))+0.8/(1+4*omw2(1,ii)));
    RR1inter=(0.1+0.3/(1+omw2(1,ii))+0.6/(1+4*omw2(1,ii)));
    RR1ampl=(f1_S1(ii)*beta_intra*tauC(1,ii)-alpha2(1,ii))/(1.0+WD);
%%%%
    RR1(1,ii)= (RR1ampl)*( RR1intra + WD*RR1inter );
%%%%
    RR2intra=(0.2/(1+omw22(1,ii))+0.8/(1+4*omw22(1,ii)));
    RR2inter=(0.1+0.3/(1+omw22(1,ii))+0.6/(1+4*omw22(1,ii)));
    RR2ampl=(f2_S2(ii)*beta_intra*tauC2(1,ii)-alpha2(1,ii))/(1.0+WD);
%%%%
    RR2(1,ii)= (RR2ampl)*( RR2intra + WD*RR2inter );
  end; %% loop ii
%%
%%      teof1=alfa2+RR1+RR2
  ymodel(1,:) = alpha2 + RR1 + RR2;
%%
end %% end R1_R2model_v1()
