function MdcevOutput(results_b, results_atrans, results_gammatrans,...
    results_scaletrans, time, sizes, names, ll, delta_method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates output for the MDCEV model
%
% Inputs:  results_b - The results from estimation
%          time - The time that optimization took to complete
%          sizes - sizes of data/parameters
%          names - The names of parameters
%          ll - The value of the log-likelihood function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all_names = strvcat('');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results Section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\n');
fprintf('**************************************************************\n')
fprintf('*********************Iterations Completed*********************\n')
fprintf('**************************************************************\n\n')

%Model information
if sizes.model_type==1
    fprintf('Model type: Gamma profile \n');
elseif sizes.model_type==2
    fprintf('Model type: Alpha profile \n');
elseif sizes.model_type==3
    fprintf('Model type: Hybrid gamma profile \n');  
end
% Scale information
if sizes.scale_type == 0
    fprintf('All scale parameters set to unity \n');  
elseif sizes.scale_type == 1
    fprintf('One scale parameter estimated \n');    
elseif sizes.scale_type == 2
    fprintf('Two scale parameters estimated \n'); 
end
fprintf('Number of observations = %4.4f \n', sizes.nobs);
fprintf('Number of estimated parameters = %4.4f \n', sizes.k);
fprintf('Total time to complete estimation was %4.4f seconds.\n', time);
fprintf('The log-likelihood value is:  %4.4f\n', ll);
fprintf('BIC =  %4.4f\n', -2*ll+log(sizes.nobs)*sizes.k);
fprintf('AIC = %4.4f\n', -2*ll+2*sizes.k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Estimates Section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n');
fprintf('**************************************************************\n')
fprintf('*********************Parameter Estimates**********************\n')
fprintf('**************************************************************\n\n')

% Index Counters
pos = 1;
psi_pos = 1;
alpha_pos = 1; 
gamma_pos = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Estimates for psi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

in.cnames = strvcat('estimate','std. err.', 't-statistic');
in.rnames = strvcat('variables', names.psi);
in.fmt    = '%6.4f';
disp('');
fprintf('RESULTS FOR PSI PARAMETERS, PSI \n');
fprintf('************************************************************\n');
mprint(results_b(pos:pos+sizes.nvar_psi-1, :),in);
disp('');

% Update Counters
pos = pos+sizes.nvar_psi;

if delta_method==0
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Parameter Estimates for alpha (raw)
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     in.cnames = strvcat('estimate','std. err.', 't-statistic');
     in.rnames = strvcat('variables', names.alpha);
     in.fmt    = '%6.4f';
     disp('');
     fprintf('RESULTS FOR SATIATION PARAMETERS, ALPHA (UNTRANSFORMED: ALPHA=1-EXP(EST)) \n');
     fprintf('************************************************************\n');
     mprint(results_b(pos:pos+sizes.ngoods, :),in);
     disp('');

     % Update Counters
     pos = pos+sizes.ngoods+1;

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Parameter Estimates for gamma
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         in.cnames = strvcat('estimate','std. err.', 't-statistic');
         in.rnames = strvcat('variables', names.gamma);
         in.fmt    = '%6.4f';
         disp('');
         fprintf('RESULTS FOR TRANSLATION PARAMETERS, GAMMA \n');
         fprintf('************************************************************\n');
     if sizes.model_type==2 
            mprint(results_b(pos:pos+sizes.ngoods-1, :),in);
            pos = pos+sizes.ngoods; % Update Counters
     else      
        mprint(results_b(pos:pos+sizes.nvar_g-1, :),in);
        pos = pos+sizes.nvar_g; % Update Counters       
     end
         disp('');

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Parameter Estimates for scale
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     in.cnames = strvcat('estimate','std. err.', 't-statistic');
     in.rnames = strvcat('variables', names.scale);
     in.fmt    = '%6.4f';
     disp('');
     fprintf('RESULTS FOR SCALE PARAMETERS (UNTRANSFORMED: SCALE = EXP(EST)) \n');
     fprintf('************************************************************\n');
     mprint(results_b(pos:end, :),in);
     disp('');

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %DISPLAY TRANSFORMED ALPHA AND SCALE PARAMETERS
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif delta_method==1

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Parameter Estimates for alpha (transformed)
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     in.cnames = strvcat('estimate','std. err.', 't-statistic');
     in.rnames = strvcat('variables', names.alpha);
     in.fmt    = '%6.4f';
     disp('');
     fprintf('RESULTS FOR SATIATION PARAMETERS, ALPHA (TRANSFORMED: ALPHA=1-EXP(EST)) \n');
     fprintf('************************************************************\n');
     mprint(results_atrans,in);
     disp('');
 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Parameter Estimates for gamma (transformed)
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     in.cnames = strvcat('estimate','std. err.', 't-statistic');
     in.rnames = strvcat('variables', names.gamma);
     in.fmt    = '%6.4f';
     disp('');
     fprintf('RESULTS FOR TRANSLATION PARAMETERS, GAMMA (TRANSFORMED: GAMMA=EXP(EST)) \n');
     fprintf('************************************************************\n');
     mprint(results_gammatrans,in);
     disp('');
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Parameter Estimates for scale
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     in.cnames = strvcat('estimate','std. err.', 't-statistic');
     in.rnames = strvcat('variables', names.scale);
     in.fmt    = '%6.4f';
     disp('');
     fprintf('RESULTS FOR SCALE PARAMETERS (TRANSFORMED: SCALE = EXP(EST)) \n');
     fprintf('************************************************************\n');
     mprint(results_scaletrans,in);
     disp('');
      
end