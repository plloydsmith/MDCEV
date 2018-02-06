function util = MdcevUtil(demand, inc, price, alpha, gamma, psi, sizes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates direct utility using the MDCEV utility function
%
% Inputs:   demand - The non-numeraire demands
%           inc - income for each individual
%           price - individual price for each good
%           alpha - The parameter estimates for alpha
%           gamma - The parameter estimates for theta
%           psi - The parameter estimates for psi
%           sizes - A structure containing information on model type
%
% Outputs:  util - The MDCEV utility
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sizes.model_type == 1
    util = (1 ./ alpha(:,1)) .* psi(:,1) .*...
        (inc - sum(price(:,2:end) .* demand(:,2:end),2)) .^ alpha(:,1)+ ... % numeraire
        sum((psi(:,2:end) .* gamma(:,2:end)) .*... 
        log(demand(:,2:end) ./ gamma(:,2:end) + 1),2);

elseif sizes.model_type ~= 1
    util = (psi(:,1) ./ alpha(:,1)) .*...
        (inc - sum(price(:,2:end) .* demand(:,2:end),2)) .^ alpha(:,1)+ ... % numeraire
        sum((psi(:,2:end) .* gamma(:,2:end) ./ alpha(:,2:end)) .*...
        ((demand(:,2:end) ./ gamma(:,2:end) + 1) .^ alpha(:,2:end) - 1),2);
end

