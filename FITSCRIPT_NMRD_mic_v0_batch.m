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
N_MCMC   = 4000;
N_sub = 100;
MCSTART=0
NBATCH=5;
beta=[0.0];
%%beta=[0.0041,0.0164,0.1];
%%beta=[1e-6,4e-6,1.6000e-05,6.4000e-05,2.5600e-04,0.0010,0.0041,0.0164,0.1,0.3,1.0];
load beta2.mat;
%%beta=[1e-6,4e-6,1.6000e-05,6.4000e-05,2.5600e-04,0.0010,0.0041,0.0164,0.1,0.3,1.0,0.02,0.12,0.22,0.32,0.06,0.15,0.45,0.6,0.8,0.01,0.035];
beta=beta2;
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
vars.seed=151;
vars.ydata=ydata;
filen='MCMC_v0c_beta';
const.startfile='MCMC_v0c_B/MCMC_v0c_beta_27.mat';
if(MCSTART==0)
   tmp=load(const.startfile);
  mcmcout=tmp.mcmcout;
  disp('starts with conf from file:');
  vars.mcmcout=mcmcout;
end;
if(NBATCH==0)
  vars.iimc=1;
  [mcmcout]  = MCMCdriver_v1(vars);
  tmpname=strcat(filen,num2str(NBATCH),'.mat');
  save(tmpname,'mcmcout');
else
%% multiple batch-jobs remotely:
   c = parcluster; % handle to the cluster

   c.AdditionalProperties.AccountName = 'nhakanss';
   t0 = tic;
   inpWallTime='04:29:00'
   disp('format of submitted jobs:         >> j_job{ii} = batch(c,@MCMCdriver_v1,1,{vars},''CurrentFolder'',RemoteDir);');
   disp('after all jobs are submitted ...  >> j_jobs{ii}.wait; %% (until finished)');
   disp('jobs are retrived with:           >> BDout = j_jobs{ii}.fetchOutputs{:}; %% ... (until no more jobs on puhti)');
   RemoteDir='/users/nhakanss/matlab_scratch/'
   req_mem=5000;
   ll = length(num2str(req_mem));
   mem_allocation = repmat(char(0),1,ll+1);
   mem_allocation(1,1:ll) = num2str(req_mem);
   mem_allocation(1,ll+1) = 'M';
   fprintf('requesting memory %s\n',mem_allocation);
   fprintf('requesting walltime %s\n',inpWallTime);
   c.AdditionalProperties.WallTime = inpWallTime;
   c.AdditionalProperties.MemUsage = mem_allocation;
   c.AdditionalProperties.QueueName = 'small';
   c.AdditionalProperties.AccountName = 'project_2002241'
   c.saveProfile;
   compl = zeros(1,NBATCH); % flag for completed job
   j_jobs= cell(1,NBATCH);  % handle of the jobs
   % submit all jobs:
   Nsub = 0;     
   for ii=1:NBATCH,
     vars.iimc=ii;
     if(compl(ii)==0 )
        if(MCSTART==0)
          %%tmpname=filen;
          %%tmpname=strcat(filen,'_',num2str(ii),'.mat');
          %%mcmcout=load(const.startfile);
          vars.mcmcout=mcmcout;
        end;
        %% set seed
        vars.seed=vars.seed+ii;
        fprintf('submitting MCMC job# %d , with seed: %d\n',ii,vars.seed);
        %%vars.mcmcout
        %% submit remotely:
        j_jobs{ii} = batch(c,@MCMCdriver_v1,1,{vars},'CurrentFolder',RemoteDir);
        Nsub = Nsub + 1;
        test_State = j_jobs{ii}.State;
        MM=0;
        while((test_State(1:3) == 'que') & (MM<5))
          pause(20); %% Wait 15s for batch queue to get organized before next job
          test_State = j_jobs{ii}.State;
          MM= MM + 1;
        end;
      end; %% if-compl(ii)
   end; %% for-ii
   fprintf('Submitted %d jobs to puhti|\n',Nsub);
%%if(NBATCH>0) 
   check_time = 2;
   while( sum( abs(compl) )<NBATCH )
     pause(check_time);
     for ii=1:NBATCH,
         if(~compl(ii)) % if job ii not registered as completed then
            %% check if completed:
            if( sum(j_jobs{ii}.State(1:3)=='fin')==3  )
              j_jobs{ii}.wait;
              % collect the result
              mcmcout = j_jobs{ii}.fetchOutputs{:};
              tmpname=strcat(filen,'_',num2str(ii),'.mat');
              save(tmpname,'mcmcout');
              fprintf('job: %d completed:-) , and saved!',ii);
              if(mod(ii,5)==0)
                fprintf('\n');
              end;
              compl(ii) = 1; % note as completed
              %% clear job ii
              j_jobs{ii}.delete;
              clear mcmcout;
            else %% job is not finished
              %% check if failed
              if( sum(j_jobs{ii}.State(1:3)=='fai')==3 )
                % stop the while loop, need debugging
                fprintf('job %d reported as failed:-( \n',ii);
                fprintf(' %s',j_jobs{ii}.Parent.getDebugLog(j_jobs{ii}.Tasks(1)));
                compl(ii) = -1;
                j_jobs{ii}.delete;
              end; %% if-fail
            end; %% if-check-completed
         end; %% if-~compl
     end; %% for-ii
   end; %% while-<NBATCH    
 end; %% if-NBATCH
return;
