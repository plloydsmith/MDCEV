%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maximum likelihood estimation of Chandra Bhat's               %
% MCDEV specification with a single essential outside good      %
% Version 1.0 - March 2017                                      %
%                                                               %
% Coded by:                                                     %
% Joshua Abbott                                                 %
% School of Sustainability                                      %
% Arizona State University                                      %
%                                                               %
% Patrick Lloyd-Smith                                           %
% Agricultural and Resource Economics                           %
% University of Saskatechewan                                   %
%                                                               %
% Elements of code by H. Allen Klaiber and Chandra Bhat are     %
% gratefully acknowledged                                       %                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The direct utility specification estimated here is of the following form:
%
% U(X)=(1/alpha1)*exp(eps1)*(y-p'x)^alpha1 + SUM((gammaK/alphaK)*PsiK*...
%{((xK/gammaK)+1)^alphaK-1})
%
% alphaK=1-exp(aK) s.t. alphaK<1
% gammaK=exp(lgammaK) s.t. gammaK>0
%
% Allow for the alphas and gammas to be heterogeneous across goods but not
% across individuals.  Parameterize PsiK as an exponential function of
% observable characteristics of the site.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT SECTION - Provide a file to save the results, a file for
% printing the output, and a title for the output.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
output_file  = 'Results\mdcev_output.out';
name_results_mat = 'Results\mdcev_results.mat';
title        = 'MDCEV model for Simulated Data';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA SECTION - The user must input or construct the following matrices, each of which
%                has rows equal to nobs (number of people in a cross section data set)
%
%   Matrices:
%       ngoods: Number of non-numeraire goods.
%       nobs:   Number of observations
%       price:  Matrix holding prices for all the goods.  Has dimensions
%               nobs x ngoods.
%       quant:  Matrix holding quantities consumed for all non-numeraire goods.  Has
%               dimensions nobs x ngoods.  Should be 0 if a censored good.
%       inc:    Vector holding annual income for all people included in the
%               model.  Has dimension nobs x 1.
%       C:      Matrix indicating which goods are available to which
%               people.  Has dimensions nobs x ngoods.  C is a matrix of 1's
%               when all goods are available to all people.
%       weight: Vector holding sample weights.  Has dimension nobs x 1.
%       het_scale: Vector holding indicators for scale.  Has dimension nobs x 1.
%
%       In addition there should be other site or individual
%       characteristics.  Site-specific variables should be in separate nobs x ngoods
%       matrices. Individual specific variables should be in nobs x 1 vectors.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set root path (change this to work for your machine and all other paths should work)
%Place this root along with all sub directories on the matlab path
cd T:\Code\MDCEV
addpath(genpath('T:\Code\MDCEV'))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Data\simulated_data.mat

% Parameter values for simulated data
% nobs = 500
% ngoods = 10
% model_type 3
% beta = [-5; 0.5; 0.5; 1.5; 3];
% gamma = 2
% alpha = 0.6

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Choose model specification =1, 2 or 3, see below
model_type=3;

% Choose scale specification
% scale_type = 0: All scales fixed
% scale_type = 1: One scale parameter estimated
% scale_type = 2: Two scales estimated
scale_type=1;

% Choose whether to restrict Alpha between 0 and 1
%alpha_rest = 0: alpha only restricted to be less than 1
%alpha_rest = 1: alpha restricted to be between 0 and 1 using logit formula
alpha_rest=0;

%Choose whether to display transformed versions of alpha and scale
%parameters w/delta method std. errors and z-statistics
delta_method=1;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLE/DATA OUTPUT NAMES SECTION - enter names for the variables
% See Bhat (2008) for descriptions of these functional forms
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psi_names = char('constant', 'attribute1',...
    'attribute2', 'socio1', 'socio2');

if scale_type == 0 || scale_type == 1
      scale_names = char('scale');
elseif scale_type == 2
    scale_names = char('scale1', 'scale2');
end

gamma_names = char('gamma2', 'gamma3', 'gamma4', 'gamma5', 'gamma6',...
    'gamma7', 'gamma8', 'gamma9', 'gamma10', 'gamma11');
alpha_names = char('alpha1', 'alpha2', 'alpha3', 'alpha4', 'alpha5',...
    'alpha6', 'alpha7', 'alpha8', 'alpha9', 'alpha10', 'alpha11');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create useful variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id = (1:1:size(quant,1))';
[nobs, ngoods] = size(quant);
C      = ones(nobs,ngoods); %Modify if censoring present
weight    = ones(nobs,1); %Modify if weights are used
het_scale    = ones(nobs,1); %Modify if hetergenous scale are used

%Create price and expenditures for the numeraire (Hicksian composite good)
price_num = ones(nobs,1);
exp_num = inc - sum(price.*quant,2);

%Create useful columns of zeros and ones
zero = zeros(nobs,1);
one = ones(nobs,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DESCRIBE MODEL SPECIFICATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SPECIFY THE VARIABLES IN EACH PSI PARAMETER
%ADD AS MANY LINES AS CHOICES
%NUMBER OF COLUMNS EQUALS THE NUMBER OF UNIQUE PARAMETERS IN THE PSI
%PARAMETERS. EACH PSIK MUST BE OF SAME DIMENSION
%DUMMY VARIABLES ARE INDICATED BY "ONE"

psi2 = [one, one, one, psi_socio];
psi3 = [one, one, zero, psi_socio];
psi4 = [one, one, one, psi_socio];
psi5 = [one, one, zero, psi_socio];
psi6 = [one, one, one, psi_socio];
psi7 = [one, zero, zero, psi_socio];
psi8 = [one, zero, one, psi_socio];
psi9 = [one, zero, zero, psi_socio];
psi10 = [one, zero, one, psi_socio];
psi11 = [one, zero, zero, psi_socio];

psi = [psi2; psi3; psi4; psi5; psi6; psi7; psi8; psi9; psi10; psi11];
nvar_psi = size(psi,2);     %Number of unique variables in psi
clear psi2 psi3 psi4 psi5 psi6 psi7 psi8 psi9 psi10 psi11;  %Keeps memory free

% APLHA and GAMMA (Don't need to set)
if model_type==1 || model_type==3
    a = ones(nobs * (ngoods+1),1);
    g = -1*eye(ngoods);
    g = repmat(g,nobs,1);
    g = -1*sortrows(g);
elseif model_type==2
    a = -1*eye(ngoods+1);
    a = repmat(a,nobs,1);
    a = -1*sortrows(a);   
    g=[];
end

nvar_a = size(a,2);
nvar_g = size(g,2);         %Number of variables in gamma

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STARTING VALUES SECTION - Enter a vector beta_start that holds starting
% values for the model.
%
% The sequence of starting values MUST BE:
%          1.  psivars
%          2.  alphas
%          3.  gammas
%          4.  scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Favored starting points
psi_start = zeros(nvar_psi,1);
a_start = zeros(nvar_a,1); %Large negative numbers are good
g_start = ones(nvar_g,1);

if scale_type == 0
    scale1_start = []; %Scale is exp(scale_start)
    scale2_start = []; 
elseif scale_type == 1
    scale1_start = 0; 
    scale2_start = []; 
elseif scale_type == 2
    scale1_start = 0; 
    scale2_start = 0;
end

%Random starting points (to test robustness)
%psi_start = rand(nvar_psi,1);
%a_start = -1*rand(nvar_a,1); %Large negative numbers are good
%g_start = rand(nvar_g,1);
%scale_start = rand(1,1); %Scale is exp(scale_start)

beta_start = [psi_start; a_start; g_start; scale1_start; scale2_start];

nvar_scale = size(beta_start,1)- nvar_psi - nvar_a - nvar_g ;%Number of variables in scale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIMIZATION SECTION - OPTIONAL
%
% This section allows the user to provide options for maximum likelihood.
% For more information, see the MATLAB help index for fminunc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt = optimset('GradObj','off',...
                           'Hessian','off',...
                           'LargeScale','off',...
                           'Display','iter',...
                           'MaxIter', 100000,...
                           'MaxFunEvals', 1000000, ...
                           'TolFun', 1.000e-006,...
                           'TolX', 1.000e-006,...
                           'FunValCheck','off',...
                           'DerivativeCheck','off',...
                           'Diagnostics','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT CHANGE ANYTHING BELOW THIS POINT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checks on Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if range(het_scale) == 0
    if scale_type == 2
        fprintf('No variation in het_scale indicator. Need to use scale_type = 0 or scale_type = 1. \n');
        return
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diary off;
if exist(output_file, 'file') ~= 0
    delete(output_file);
end
diary(output_file);

fprintf(1, title);
fprintf(1, '\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store the data in a structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = struct('psi', psi, 'price', price, 'quant', quant,...
   'inc', inc, 'C', C, 'weight', weight, 'price_num', price_num, 'exp_num', exp_num,...
   'a', a, 'g', g, 'het_scale', het_scale);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store variable names in a structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
names = struct('scale', scale_names, 'psi', psi_names, 'gamma', gamma_names,...
    'alpha', alpha_names);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine sizes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[k, junk] = size(beta_start); % Determine number of params

% Define a struct to hold sizes and model type
sizes = struct('nvar_psi', nvar_psi,'nvar_g', nvar_g, 'nvar_a', ...
    nvar_a, 'nvar_scale', nvar_scale, 'ngoods', ngoods, 'nobs', nobs,...
    'model_type', model_type, 'k', k, 'scale_type', scale_type,...
    'alpha_rest', alpha_rest);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Index structure for coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ndx = MdcevNdx(sizes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Maximum Likelihood estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time = clock;
[b_kt,fval,exitflag,output,grad,hessian] = fminunc(@(b) MdcevLike(b, ...
   data, sizes, ndx), beta_start, opt);
time = etime(clock,time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate log likelihood values
ll = -1*MdcevLike(b_kt, data, sizes, ndx);

% Caculate cov_mat, standard errors, and t-statistics
cov_mat = inv(hessian);
std_err = sqrt(diag(cov_mat));
zstat = b_kt./std_err;

% Save results in results matrix
results_b = [b_kt std_err zstat];

% Reconfigure results_b matrix
results_psi = results_b(1:ndx.psi(end),:);

% Add back in other alphas for output/export
results_alpha = results_b(ndx.alpha(1):ndx.alpha(end),:);

if model_type ==1
    results_alpha = [results_alpha; repmat([0, 0, 0], ngoods, 1)];
elseif model_type ==3
    results_alpha = repmat(results_alpha, ngoods+1,1);
end

% Add back in other gammas for output/export
if model_type==2
    results_gamma = repmat([0, 0, 0], ngoods,1);
elseif model_type~=2
    results_gamma = results_b(ndx.gamma(1):ndx.gamma(end),:);
end

% Add in other scale parameter if fixed
results_scale = results_b(end-nvar_scale+1:end,:); %Two scales

if scale_type == 0
    results_scale = [results_scale; 0, 0, 0];
end

results_b = [results_psi; results_alpha; results_gamma; results_scale];

% Transform estimates for alpha and scale (treating parameters as
% coefficients on mutually exclusive dummy variables)
astart=sizes.nvar_psi+1;
alpha_results = results_b(astart:astart+sizes.ngoods,:);

if sizes.alpha_rest == 0
    a1 = 1-exp(alpha_results(:,1)); %Transform estimates to alphaKs
    a2 = exp(alpha_results(:,1)).*alpha_results(:,2); %Delta method std. error
elseif sizes.alpha_rest == 1
    a1 = exp(alpha_results(:,1))/(1+exp(alpha_results(:,1))) ; %Transform estimates to alphaKs
    a2 = exp(alpha_results(:,1)).*alpha_results(:,2); %Delta method std. error
end

a3 = a1(:,1)./a2;
results_atrans = [a1(:,1) a2 a3];
clear alpha_results a1 a2 a3;

s1 = exp(results_scale(:,1)); %Transform estimates to scale parameters
s2 = exp(results_scale(:,1)).*results_scale(:,2); %Delta method std. error
s3 = s1./s2;

results_scaletrans = [s1 s2 s3];
clear results_scale s1 s2 s3;

g1 = exp(results_gamma(:,1)); %Transform estimates to scale parameters
g2 = exp(results_gamma(:,1)).*results_gamma(:,2); %Delta method std. error
g3 = g1./g2;

results_gammatrans = [g1 g2 g3];
clear results_gamma g1 g2 g3;

% Print the output
MdcevOutput(results_b, results_atrans, results_gammatrans,...
    results_scaletrans, time, sizes, names, ll, delta_method)
diary off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta=results_b(ndx.psi(1):ndx.psi(end), 1);
est=struct('alpha', results_atrans(:,1), 'scale', results_scaletrans(:,1),...
    'beta', beta, 'gamma', results_gammatrans(:,1));

save(name_results_mat, 'data', 'names', 'ndx', 'sizes', 'est', 'id');
