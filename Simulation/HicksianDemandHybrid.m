function [hdemand] = HicksianDemandHybrid(util, price, MUzero, alpha, gamma, sizes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates hicksian demands for a given dataset and set of error
% terms assuming Bhat's (2008) MCDEV 'gamma profile' specification model.  
% See Lloyd-Smith (2017) for explanation of steps.
% 
%   Outputs: 
%       hdemand - nobs X ngoods array of Hicksian demands  
%
%   Inputs:
%       util - value of utility in basecase
%       price - individual price for each good
%       MUzero - price normalized MU at zero consumption
%       alpha, gamma -  at altered bases
%       sizes - structures defined as in EstimateMDCEV.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%j=56;alpha = alpha_b;gamma = gamma_b;MUzero = MUzero_b;price = price_p;inc = inc_b;psi = psi_b0;util = util;
hdemand=zeros(sizes.nobs,sizes.ngoods+1); 

%Cycle over each observation   
for j=1:sizes.nobs
    
    M=1; %Indicator of which ordered alteratives (<=M) are being considered  
	exit=0;
	A=[MUzero(j,:)', gamma(j,:)', price(j,:)', alpha(j,:)', (1:sizes.ngoods+1)']; %Make matrix of MU at zero, gamma, price and good indicator
	B=sortrows(A(2:end,:),-1); %Sort matrix in descending order of non-numeraire utilities
	C=[A(1,:);B]; %Add essential numeraire at the top of the rank ordering
	e=C(:,1).*C(:,3); %%%%%%%%%% obtain psi term
	a=C(:,2).*e; %%%%%%%%%% Want psi not price term
	b=C(:,1).^(-C(:,4)./(C(:,4)-1)); %%%%%%%%%% want price/psi so take negative of exponent
	c=a.*b;

    while exit==0
        %Calculate 1/lambda for a given M
        lambda_num=C(1,4)*util(j)+sum(a(1:M,:),1)-e(1,:); %%%%%%%%%% utility not expenditure and subtract numeriare psi
        lambda_den=sum(c(1:M,:),1);
        lambda=(lambda_num/lambda_den).^(-(C(1,4)-1)/C(1,4)); %%%%%%%%%% new exponent term
        lambda1 = 1/lambda;    %%%%Compare 1/lambda to MU
            
        %Compare 1/lambda to baseline utility of the next lowest alternative
        %(M+1). If lambda exceeds this value then all lower-valued
        %alternatives have zero demand.
        if (lambda1 > C(min(M+1,sizes.ngoods+1),1) || M==sizes.ngoods+1)
            %Compute demands (using eq. 12 in Pinjari and Bhat)
            d=[0;ones(sizes.ngoods,1)];
            % Need 'real' command as sometimes exponentiate negative demand
            X=real( ((1./(lambda.*C(:,1))).^(1./(C(1,4)-1))-d).*C(:,2)); 
            X=[X(1:M,:);zeros(sizes.ngoods+1-M,1)];
            exit=1;
        elseif M < sizes.ngoods+1
            M=M+1;           
        else                     
        end
    end
    
    %This code puts the choices back in their original order and
    %exports demands
    e=[X,C(:,5)];
    f=sortrows(e,2);
    hdemand(j,:)= f(:,1)';        
        
end
    