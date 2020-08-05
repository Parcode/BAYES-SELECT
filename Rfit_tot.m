%% call: [err] = Rfit_tot(x,ydata,const,yerr);
%%
%% Compute chi2 from experimental data "ydata" and given model "R1_R2_model".
%% Note, the actual model function is passed as argument "R1_R2model" with additional parameters passed with "const",
%% making Rfit_tot a general format for most problems.
%%
%% Input: x: array of fitting parameters
%%        ydata: experimental relaxation rates (DeltaR = 1/T2 - 1/T1 )
%%        const: fixed constants in the model (translational diffusion ...)
%%               and other constants like h-bar , gamma, etc. needed in R1_R2_model.
%%        yerr: standard deviation of ydata
%% Output: err = sum_i (ydata_i - model_i)/yerr_i )
%% Note the call of R1_R2model has the arguments: [ymodel]=R1_R2model(ydata,const,jj); 
%% where jj is the experimental point number jj.
%%
function [errtot,err,ymodel] = Rfit_tot(x,ydata,const,yerr);
%%
  err = zeros(const.NN,1);
  errtot=0.0;
  switch const.Nfunc
    case 0
      R1_R2model=@R1_R2model_v0c;
    case 6
      R1_R2model=@R1_R2model_v6;
    case 10
      R1_R2model=@R1_R2model_v10c;
    case 11
      R1_R2model=@R1_R2model_v11c;
    otherwise
      return;
  end;
%%  R1_R2model=const.R1_R2model;
  ymodel = zeros(const.NN,const.NNtypes); 
%%
  for jj=1:const.NN,
    [ymodel(jj,1:const.NNtypes)] = R1_R2model(x,ydata,const,jj);
    for kk=1:const.NNtypes,
      err(jj,1) = err(jj,1) + (ydata(jj,kk)-ymodel(jj,kk) )^2/yerr(jj,kk)^2;
    end;
  end;
  errtot = sum(err(:));
%% errtot = 0.5*errtot/(Nobs-Nfit);
  errtot = (const.beta)*0.5*errtot;
end

