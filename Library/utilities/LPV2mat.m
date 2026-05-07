%% LPV2mat
% ------------------------------------------------------------------------
% Construction of the map M: R -> R^(mxn)
% 
% Pablo R Baldivieso Monasterios
% 21/09/2022
% The University of Sheffield
% ------------------------------------------------------------------------

function A = LPV2mat(M,rho)

rho0 = M.bk3.first;
drho = M.bk3.spacing;
Nrho = M.bk3.points;

n = M.bk1.points;
m = M.bk2.points;

rho_points = ((0:Nrho-1)*drho + rho0)';

A = zeros(n,m);

for i = 1:n 
    for j = 1:m 
        A(i,j) = interp1(rho_points,squeeze(M.dat(i,j,:)),rho);
    end
end