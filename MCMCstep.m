% MCMC with Metropolis-Hasting algorithm. parameters are bounded by pmin and pmax.
% CALL: [Xout, acc, newXbest,errbest,err] = MCMCstep(Xout, yerr.^2, const,h,xmin,xmax, Xbest,ydata,errbest,@Rfit_tot_D);
%
% INPUT: X0: 
% 
function [Xout, acc, newXbest,errbest,err] = MCMCstep(X0, STD, const,h,pmin,pmax, Xbest,ydata,errbest,mymodel)
% number of model configurations (parameters):
Nparm = length(X0);
yerr = STD; %% Rfit_tot_C takes the standard deviation as input
%
%% uniform random number for MC-step ( mean zero, var(xi)=1.) 
xi = sqrt(3.0).*(2.*rand(1,Nparm) - 1.0);
%% individual configuration steplength, xi scaled by h:
%%xi
%%h
dxi = (h.*xi);
%
% Trial step Xtr
Xtr = X0;
% save the previos step
X0tmp = X0;
%%function [errtot,err,ymodel] = Rfit_tot(x,ydata,const,yerr,R1_R2model);
[errtmp] = mymodel(X0tmp,ydata,const,yerr);
% book-keep what steps was accepted (acc)
acc = zeros(1,Nparm);
% book-keep the currently best configuration (Xbest)
% note however, that in a baysian sence the most likely configuration set 
% is computed from the probability distribution
newXbest = Xbest;
%
for ii = 1:Nparm,
%%for ii = indNparm,
% make the step in configuration space
  Xtr(1,ii) = X0tmp(1,ii) + dxi(1,ii);
%% if bounded domain, reflections at boundary:
  if(Xtr(1,ii) < pmin(1,ii) )
    Xtr(1,ii) = 2*pmin(1,ii) - Xtr(1,ii);
  end; % "if min"
  if (Xtr(1,ii) > pmax(1,ii) )
    Xtr(1,ii) = 2*pmax(1,ii) - Xtr(1,ii);
  end; % "if max"
%% Second round (in case large dxi) 
  if(Xtr(1,ii) < pmin(1,ii) )
    Xtr(1,ii) = 2*pmin(1,ii) - Xtr(1,ii);
  end; % "if min"
  if (Xtr(1,ii) > pmax(1,ii) )
    Xtr(1,ii) = 2*pmax(1,ii) - Xtr(1,ii);
  end; % "if max"
%% third round (in case large dxi) 
  if(Xtr(1,ii) < pmin(1,ii) )
    Xtr(1,ii) = 2*pmin(1,ii) - Xtr(1,ii);
  end; % "if min"
  if (Xtr(1,ii) > pmax(1,ii) )
    Xtr(1,ii) = 2*pmax(1,ii) - Xtr(1,ii);
  end; % "if max"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute mean square error for trial move  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [err] = mymodel(Xtr,ydata,const,yerr);
  Lfrac = exp( -(err - errtmp) );
  acc(1,ii) = min([1.0,Lfrac]);
  u = rand(1,1); % uniform random number [0,1]
  if(u < acc(1,ii) ) %% then accept the configuration move
    %% keep the updated Xtr(ii,1)
    X0tmp(1,ii) = Xtr(1,ii);
    acc(1,ii) = 1;
  else % keep the old configuration
    Xtr(1,ii) = X0tmp(1,ii); % move is rejected)
    acc(1,ii) = 0;
  end; % "if u"
  if(err<errbest)
    fprintf(' new lowest err :-) %5.3e old: %5.3e\n',err,errbest);
    errbest = err;
    newXbest=Xtr;
  else
    newXbest = Xbest;
  end;
end; % for - Nparm
Xout = Xtr;
end % MCMCstep

