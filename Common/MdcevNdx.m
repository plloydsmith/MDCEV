function ndx = MdcevNdx(sizes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns a struct of indices for psi, alpha, gamma and scale
%
% Returns:  ndx.psi, ndx.gamma, ndx.alpha, ndx.scale
%
% Inputs:  sizes - A structure of sizes of data and model type
%
% Outputs:  ndx - A struct listing indices of where in the parameter vector
%                 variables start and end 
%
% Written by Josh Abbott 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create index of good specific effects, psi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Find starting and ending points
b_psi = 1;
e_psi = (b_psi-1) + sizes.nvar_psi;

% create index
ndx_psi = [b_psi e_psi];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create index of alphas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find starting and ending points
b_alpha = e_psi+1;
e_alpha = (b_alpha-1) + sizes.nvar_a;

% Create index
if b_alpha==e_alpha %Allows for possibility of single alpha (gamma profile, hybrid profile)
    ndx_alpha=[b_alpha];
else
ndx_alpha = [b_alpha e_alpha];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create gamma index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sizes.model_type~=2   
    %i.e. not an alpha profile model%
    b_gamma = e_alpha+1; 
    e_gamma = (b_gamma-1) + sizes.nvar_g;
    
    %Create index
    ndx_gamma = [b_gamma e_gamma];
else
    ndx_gamma = [];
    e_gamma=e_alpha; %So scale code works
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create scale index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b_scale = e_gamma+1; 
e_scale = (b_scale-1) + sizes.nvar_scale;
ndx_scale = [b_scale e_scale];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create struct to return indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ndx = struct('psi', ndx_psi, 'alpha', ndx_alpha, 'gamma', ndx_gamma, 'scale', ndx_scale);