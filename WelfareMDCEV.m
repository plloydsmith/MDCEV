%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WelfareMDCEV.m --
%Program to Compute Hicksian/Marshallian Demands and Welfare for MDCEV
%model
%
% Coded by:                                                     %
% Joshua Abbott                                                 %
% School of Sustainability                                      %
% Arizona State University                                      %
%                                                               %
% Patrick Lloyd-Smith                                           %
% Agricultural and Resource Economics                           %
% University of Saskatechewan                                   %
%                                                               %                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

%Load structures from model estimation. These include the estimation
%sample.  
load Results\mdcev_results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error Draws and Simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set number of error draws 
L=25; % Number of error draws

rng(2018);

% Choose whether to draw errors conditional on actual demand or not
% 1 = Yes, 0 = No
cond_error = 1;

% Use general welfare algorthim
% 1 = Yes, 0 = No
% Note: Only model_type=3 can be used with hybrid (non-general) routine
algo_gen = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

null=zeros(sizes.nobs,1); %a handy column of zeros
one=ones(sizes.nobs,1); %a handy column of ones
mdemand=zeros(sizes.nobs,sizes.ngoods+1,L);
wtp1 = zeros(sizes.nobs,L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set up baseline data (_b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_b = bsxfun(@times, ones(sizes.nobs, sizes.ngoods+1), est.alpha');
price_b=[one, data.price]; %add numeraire price to price matrix (=1)
inc_b = data.inc;

if sizes.scale_type == 2
    scale_b = data.het_scale * est.scale(1) + (1-data.het_scale) * est.scale(2);       
else
    scale_b =  ones(sizes.nobs,1).*est.scale;  
end

gamma_b = bsxfun(@times, ones(sizes.nobs, sizes.ngoods), est.gamma');    
gamma_b = [one, gamma_b]; %add numeraire gamma parameter to gamma vector (=1)

psi_b=reshape(data.psi*est.beta, sizes.nobs, sizes.ngoods); 
psi_b = [null, psi_b];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Policy Change Scenario (_p)
% Change this section for specific policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Print out policy welfare information
policy_names = char('1: $10 price increase');

price_p = [ price_b(:,1), price_b(:,2:11)+10]; % Price change

psi_p = psi_b;
dataP=data;
gamma_p=gamma_b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sizes.model_type ~= 3 && algo_gen == 0
  fprintf('Cannot run Model with HicksianDemandGeneral approach. Choose algo_gen=1. \n');
  return
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulate EV1 errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error = DrawError(price_b, inc_b, psi_b, gamma_b, alpha_b, scale_b,... 
    data.quant, sizes, cond_error, L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimate Demands and Welfare
time = clock;

% Can set to parfor loop if necessary
for i_err=1:L
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Baseline Calculations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create psi_b w/ error term
    psi_b0 = exp(psi_b + error(:,:,i_err));    
    MUzero_b=psi_b0./price_b; % price normalized MU at zero

    if cond_error == 0
        % Estimate Unconditional Marshallian Demand for Baseline
        if algo_gen == 0
            d = MarshallianDemandHybrid(inc_b, price_b, MUzero_b, alpha_b, gamma_b, sizes);
        elseif algo_gen == 1
            d = MarshallianDemandGeneral(inc_b, price_b,MUzero_b, alpha_b, gamma_b, sizes);
        end
    elseif cond_error == 1 % Use Conditional Demands
        z = inc_b - sum(price_b(:,2:end) .* data.quant,2);
        d = [z data.quant];
    end
    
    mdemand(:,:,i_err)=d;      

    % Estimate Base Utility    
    util = MdcevUtil(d, inc_b, price_b, alpha_b, gamma_b, psi_b0, sizes);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Policy 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create new psi_p
    psi_p0 =exp(psi_p + error(:,:,i_err));    
    MUzero_p=psi_p0 ./price_p; % price normalized MU at zero

    if algo_gen == 0
       % Estimate Hicksian Demand
        h = HicksianDemandHybrid(util, price_p, MUzero_p, alpha_b, gamma_p, sizes);
    elseif algo_gen == 1
        h = HicksianDemandGeneral(util, price_p, MUzero_p, alpha_b, gamma_p, sizes);
    end

    % Calculate WTP (Policy 1)
    wtp1(:,i_err) = data.inc - sum(price_p .* h, 2);

end

time = etime(clock,time);
fprintf('SIMULATION TIME (SECS) =  %4.4f \n', time);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print Welfare Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wtp1_mean = mean(wtp1,2);

wtp_mean = mean(wtp1_mean,1);

% Print out policy welfare information
inw.cnames = char('Mean');
inw.rnames = char('Policy', policy_names(:,:));
inw.fmt = '%8.4f';
fprintf('Mean WTP for selected policies\n');
fprintf('************************************************************\n');
mprint(wtp_mean', inw);
