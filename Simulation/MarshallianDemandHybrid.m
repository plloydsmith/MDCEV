function [mdemand] = MarshallianDemandHybrid(inc, price, MUzero, alpha, gamma, sizes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function simulates marshallian demands for a given dataset and set of error
% terms assuming Bhat's (2008) MCDEV 'gamma profile' specification model.  
% See Pinjari and Bhat (2011) for explanation of steps.
% 
%   Outputs: 
%       demand - nobs X ngoods array of predictions    
%
%   Inputs:
%       inc - income for each individual
%       price - individual price for each good
%       MUzero - price normalized MU at zero consumption
%       alpha, gamma -  at altered bases
%       sizes - structures defined as in EstimateMDCEV.m
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mdemand=zeros(sizes.nobs,sizes.ngoods+1); 

for j=1:sizes.nobs
    
    M=1; %Indicator of which ordered alteratives (<=M) are being considered  
    exit=0;
    A=[MUzero(j,:)', gamma(j,:)', price(j,:)', alpha(j,:)', (1:sizes.ngoods+1)']; %Make matrix of MU at zero, gamma, price and good indicator
    B=sortrows(A(2:end,:),-1); %Sort matrix in descending order of non-numeraire utilities
    C=[A(1,:);B]; %Add essential numeraire at the top of the rank ordering
    a=C(:,2).*C(:,3); 
    b=C(:,1).^(1./(1-C(:,4)));
    c=a.*b;

    while exit==0
        %Calculate lambda for a given M
        lambda_num=inc(j)+sum(a(1:M,:),1)-1;
        lambda_den=sum(c(1:M,:),1);
        lambda=(lambda_num/lambda_den).^(C(1,4)-1);
                             
    %Compare lambda to baseline utility of the next lowest alternative
	%(M+1). If lambda exceeds this value then all lower-valued
	%alternatives have zero demand.
         if (lambda > C(min(M+1,sizes.ngoods+1),1) || M==sizes.ngoods+1)
            %Compute demands (using eq. 12 in Pinjari and Bhat)
            d=[0;ones(sizes.ngoods,1)];
            X=((lambda./C(:,1)).^(1./(C(:,4)-1)) -d).*C(:,2);
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
    mdemand(j,:)=f(:,1)';      
end
