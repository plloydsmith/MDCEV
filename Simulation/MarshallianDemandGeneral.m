function [mdemand] = MarshallianDemandGeneral(inc, price, MUzero, alpha, gamma, sizes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function simulates marshallian demands for a given dataset and set of error
% terms that can be used with a broad set of utility specifications.
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
    b=(1 ./ (C(:,4)-1));
    c=C(:,1) .^ b;
    d=[0;ones(sizes.ngoods,1)];
    inc_j = inc(j);

    % Set tolerance
    %tol_l = 1e-20;
    tol_e = 1e-20;
    max_loop = 999999;

    while exit==0
        %Calculate lambda equal to MUzero(M+1)
        lambda = C(M+1,1);
        % Calculate E
        E =  sum(a(1:M).*(lambda .^ b(1:M) ./ c(1:M)-d(1:M)),1);   
        
        if E >= inc_j || M+1 == sizes.ngoods+1
            M = (M+1).* (E < inc_j) + M .* (E >= inc_j);
            lambda_l = 0 .* (E < inc_j) + lambda.* (E >= inc_j);
            lambda_u = C(M,1);   
            lambda = (lambda_l+lambda_u) / 2;
            
            for i=1:max_loop              
                % Compute new E
                E =  sum(a(1:M) .*(lambda .^ b(1:M) ./ c(1:M)-d(1:M)),1);   
                % Update lambdas's
                lambda_u = ((lambda_l+lambda_u) / 2).*(E <= inc_j) + lambda_u.*(E > inc_j);
                lambda_l = ((lambda_l+lambda_u) / 2).*(E >= inc_j) + lambda_l.*(E < inc_j);
                lambda = (lambda_l+lambda_u) / 2;

                if abs(E-inc_j) < tol_e % || abs(lambda_l-lambda_u) < tol_l
                    break
                end     
            end    
             %Compute demands (using eq. 12 in Pinjari and Bhat)
            X=((lambda.*(1./C(:,1))).^ b -d).*C(:,2);
            X=[X(1:M,:);zeros(sizes.ngoods+1-M,1)];
            exit=1;
           
        elseif (E < inc_j)  && (M+1 < sizes.ngoods+1)
            M=M+1;             
        else           
        end
    end
        
    %This code puts the choices back in their original order and exports demands
    e=[X,C(:,5)];
    f=sortrows(e,2);
    mdemand(j,:)=f(:,1)'; 
end
