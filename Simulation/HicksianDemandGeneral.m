function [hdemand] = HicksianDemandGeneral(util, price, MUzero, alpha, gamma, sizes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function simulates hicksian demands for a given dataset and set of error
% terms that can be used with a broad set of utility specifications.
% See Lloyd-Smith (2017) for explanation of steps.
%
%   Outputs:
%       hdemand - nobs X ngoods array of predictions
%
%   Inputs:
%       util - value of utility in basecase
%       price - individual price for each good
%       MUzero - price normalized MU at zero consumption
%       alpha, gamma -  at altered bases
%       sizes - structures defined as in EstimateMDCEV.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hdemand=zeros(sizes.nobs,sizes.ngoods+1);

%Cycle over each observation
for j=1:sizes.nobs

    M=1; %Indicator of which ordered alteratives (<=M) are being considered
    exit=0;
    A=[MUzero(j,:)', gamma(j,:)', price(j,:)', alpha(j,:)', (1:sizes.ngoods+1)']; %Make matrix of MU at zero, gamma, price and good indicator
    B=sortrows(A(2:end,:),-1); %Sort matrix in descending order of non-numeraire utilities
    C=[A(1,:);B]; %Add essential numeraire at the top of the rank ordering
    b=C(:,1).*C(:,3); % obtain psi
	c=C(:,2).*b ./C(:,4) ; % gamma*psi/alpha
    d = (C(:,4)./(C(:,4)-1)); % (alpha/(alpha-1))
	e= (1./C(:,1)).^d; % (1/MUzero)^(alpha/(alpha-1))
    f=[0;ones(sizes.ngoods,1)];
    g=C(:,2).*C(:,3);

    util_j = util(j);

    % Set tolerance as a percentage change
    % value of lambda can vary quite a lot between model_types
    tol_l = 1e-20;
    %tol_u = 1e-10; only use lambda tolerance
    max_loop = 999999;

    while exit==0
        %Calculate 1/lambda equal to MUzero(M+1)
        lambda1 = C(M+1,1);

        % Calculate new utility
        util_new = sum((c(1:M) .*((lambda1 .^ d(1:M)) .* e(1:M)-f(1:M))),1);

        if util_new >= util_j || M+1 == sizes.ngoods+1
            M = (M+1).* (util_new < util_j) + M .* (util_new >= util_j);
            lambda_l = lambda1.* (util_new >= util_j) +  0 .* (util_new < util_j);
            lambda_u = C(M,1);
            lambda1 = (lambda_l+lambda_u) / 2;

            for i=1:max_loop
                % Calculate new utility
                util_new = sum((c(1:M) .*((lambda1 .^ d(1:M)) .* e(1:M)-f(1:M))),1);

                % Update lambdas's
                lambda_u = ((lambda_l+lambda_u) / 2).*(util_new <= util_j) + lambda_u.*(util_new > util_j);
                lambda_l = ((lambda_l+lambda_u) / 2).*(util_new >= util_j) + lambda_l.*(util_new < util_j);
                lambda1 = (lambda_l+lambda_u) / 2;
                if abs(lambda_l-lambda_u)/lambda_u < tol_l %|| abs(util_new-util_j)/util_j < tol_u
                    break
                end
            end
             %Compute demands (using eq. 12 in Pinjari and Bhat)
            X=((lambda1.*(1./C(:,1))).^ (1./(C(:,4)-1)) -f).*C(:,2);
            X=[X(1:M,:);zeros(sizes.ngoods+1-M,1)];
            exit=1;

        elseif (util_new < util_j)  && (M+1 < sizes.ngoods+1)
            M=M+1;
        else
        end
    end

    %This code puts the choices back in their original order and exports demands
    e=[X,C(:,5)];
    f=sortrows(e,2);
    hdemand(j,:)=f(:,1)';
end
