function error = DrawError(price, inc, psi, gamma, alpha, scale, quant, sizes, cond_error, L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function draws either unconditional or conditional error terms
%
% Inputs:   normal MDCEV terms...
%           cond_error - 1 = conditional error draw, 0 = unconditional
%           L - number of error draws
%
% Outputs:  error - The error terms
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Unconditional Error Draws
if cond_error == 0    
    error = rand(sizes.nobs, sizes.ngoods+1, L); %uniform(0,1) draws
    error = -log(-log(error)) .* scale; 
 
% Conditional Error Draws    
elseif cond_error == 1
    % Cacluate the demand function, g   
    z = inc - sum(price(:,2:end) .* quant, 2);
    cond_demand = [z quant];
   
    % compute vk and v1
    if sizes.model_type == 1
        vk = psi(:,2:end) - log(price(:,2:end)) - ...
        log(cond_demand(:,2:end) ./ gamma(:,2:end) + 1);

    elseif sizes.model_type ~= 1
        vk = psi(:,2:end) - log(price(:,2:end)) + ...
        (alpha(:,2:end) - 1).*log(cond_demand(:,2:end) ./ gamma(:,2:end) + 1);  
    end
        v1 = (alpha(:,1) - 1).*log(z);     
   
    % ek = v1 - vk and assume error term is zero for outside good
    e = (v1 - vk) ./ scale;
    e = [zeros(sizes.nobs,1), e]; 
    
    % replicate matrices for each set of error
    e = repmat(e,[1 1 L]) ;
    cond_demand_error = repmat(cond_demand,[1 1 L]) ;
    
    % Calculate errors
    % For unvisited alternatives, draw from truncated multivariate
    % logistic distribution
    error = (cond_demand_error > 0).*(e .* scale) + ...
    (cond_demand_error == 0).*((-log(-log(rand(sizes.nobs,sizes.ngoods+1,L)...
    .* exp(-exp(-e))))).* scale);       
end    
