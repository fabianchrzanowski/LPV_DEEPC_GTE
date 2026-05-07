% ========================================================================
% MATRIX INTERPOLATION
% 
% Function to interpolate a page matrix in the third dimension
%
% INPUTS
%   x: list of values in the thrid dimension
%   M: Matrix page object
%   xi: query points
%
% OUTPUTS
%   M: Interpolated matrix
%
% ========================================================================
% Hadriano Morales Escamilla
% The University of Sheffield
% 11/05/2022

function M = interp_matrix(x,M,xi)

    % bringing 3 dimension first, to be interpolated
    M = permute(M,[3,1,2]);
    
    % interpolation
    M = interp1(x,M,xi);

    % permutation back to original position
    M = permute(M,[2,3,1]);

end

