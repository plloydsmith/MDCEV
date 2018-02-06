function ll = MdcevLike(b, data, sizes, ndx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the log-likelihood of the Chandra Bhat Kuhn-Tucker Model
% 
% Inputs:  b:  vector of parameters
%               order:
%                1.  psivars  
%                2.  alphas
%                3.  gammas
%                4.  scale 
%
%          data:  struct containing data matrices for parameter types:
%               - data.psi:  psi variables
%               - data.price:  price data for each good
%               - data.quant:  quantities demanded for each good
%               - data.exp: expenditures on each good
%               - data.price_num:  numeraire price 
%               - data.exp_num:  expenditure on numeraire
%               - data.inc:  Income 
%               - data.C:  censoring information
%               - data.wgt:  sampling weights
%               - data.a: covariates for alpha
%               - data.g: covariates for gamma    
%               - data.het_scale: indicator vector for heterogenous scale
%
%          sizes:  struct contains sizes of variables and data
%               - sizes.nvar_psi: number of variables in psi
%               - sizes.nvar_g: number of gamma parameters
%               - sizes.nvar_a: number of alpha parameters
%               - sizes.ngoods:  number of goods (except numeraire)
%               - sizes.nobs:  number of choice occasions
%               - sizes.model_type: id of model (=1 gamma profile, =2 alpha
%                 profile, =3 quasi-LES)
%               - sizes.k: number of estimated parameters
%               - sizes.scale_type: type of scale
%               - sizes.alpha_rest: 1 indicates alpha constrained to be
%               between (0,1)
%
%          ndx: struct contains index of variables and parameters 
%               - ndx.psi: indices of psi parameters in b
%               - ndx.alpha: indices of alpha parameters in b
%               - ndx.gamma: indices of gamma parameters in b
%               - ndx.scale: indices of scale parameter(s) in b
%
%
% Outputs:  ll - The value of the log-likelihood function
%
% Author:  Joshua Abbott with thanks to Allen Klaiber for some code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
% Initialize likelihood
%%%%%%%%%%%%%%%%%%%%%%
like = zeros(sizes.nobs,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decompose b into its component pieces.
% These include:  scale parameter: scale
%                 site/individual characteristics: lpsi 
%                 translating parameters: gamma
%                 satiation parameters: alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%b = beta_start;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Parameters lpsi, alpha, gamma, scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lpsi = reshape(data.psi*b(ndx.psi(1):ndx.psi(end)), sizes.nobs, sizes.ngoods);

if sizes.alpha_rest == 0
    alpha = reshape(1-exp(data.a*b(ndx.alpha(1):ndx.alpha(end))), sizes.nobs, sizes.ngoods+1); 
elseif sizes.alpha_rest == 1     
    alpha1 = reshape(exp(data.a*b(ndx.alpha(1):ndx.alpha(end))), sizes.nobs, sizes.ngoods+1); 
    alpha2 = reshape(1+exp(data.a*b(ndx.alpha(1):ndx.alpha(end))), sizes.nobs, sizes.ngoods+1); 
    alpha = bsxfun(@rdivide, alpha1, alpha2); 
end

if sizes.scale_type==0
    scale = ones(sizes.nobs, sizes.ngoods+1);
elseif sizes.scale_type==1
    scale_1 = ones(sizes.nobs, sizes.ngoods+1).*exp(b(ndx.scale(1)));  
    scale_2 = ones(sizes.nobs, sizes.ngoods+1);
    scale = bsxfun(@times, data.het_scale, scale_1) + bsxfun(@times, (1-data.het_scale), scale_2); 
elseif sizes.scale_type==2
    scale_1 = ones(sizes.nobs, sizes.ngoods+1).*exp(b(ndx.scale(1)));  
    scale_2 = ones(sizes.nobs, sizes.ngoods+1).*exp(b(ndx.scale(end)));
    scale = bsxfun(@times, data.het_scale, scale_1) + bsxfun(@times, (1-data.het_scale), scale_2); 
end

if sizes.model_type==2 % Alpha profile - set all gammaK=1
    gamma = ones(sizes.nobs, sizes.ngoods);   
elseif sizes.model_type==1 || sizes.model_type==3% Gamma/Hybrid Profile - construct gammas
    gamma = reshape(exp(data.g*b(ndx.gamma(1):ndx.gamma(end))), sizes.nobs, sizes.ngoods);   
end
if sizes.model_type==1    %Gamma profile - set non-numeraire alphaK=0
    alpha = [alpha(:,1), zeros(sizes.nobs, sizes.ngoods)]; 
end
    
x=[data.exp_num./data.price_num, data.quant]; %Full demand matrix including essential numeraire
gamma_full = [zeros(sizes.nobs,1), gamma]; %Includes column for numeraire (gamma1=0)
p=[data.price_num, data.price]; %Full price matrix (including numeraire)
nonzero=(x~=0); %Logical (0,1) matrix indicating nonzero consumption
M=sum(nonzero,2); %Number of consumed goods (including numeraire)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate the determinant of the Jacobian of transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f=((1-alpha)./(x+gamma_full));
pf=sum(nonzero.*(p./f),2);
f(x==0)=1;
prodf=prod(f,2);
jacobian=prodf.*pf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate the likelihood for each observation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Log FOC terms for each good
v=lpsi+(alpha(:,2:end)-1).*log((data.quant./gamma)+1)-log(data.price);
v=[(alpha(:,1)-1).*log(data.exp_num./data.price_num), v]; %Add numeraire 
v=v./scale; %Divide by scale
v=exp(v); %Exponentiate

sumv=sum(v,2);
v(x==0)=1;
prodv=prod(v,2);
like=(prodv./(sumv.^M));
like=like.*factorial(M-1)./(scale(:,1).^(M-1)); %Must alter this for heterogenous scale
like=like.*jacobian;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Take the (negative) log of the likelihood and return that value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ll=-sum(log(like));

end
