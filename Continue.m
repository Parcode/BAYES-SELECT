   while( sum( abs(compl) )<NBATCH ) 
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

