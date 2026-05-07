%% LPV2mat
% ------------------------------------------------------------------------
% Construction of the map x0: R -> R^(n)
% 
% Pablo R Baldivieso Monasterios
% 21/09/2022
% The University of Sheffield
% ------------------------------------------------------------------------

function x0 = off2mat(x,rho)

n = x.bk1.points;

rho0 = x.bk2.first;
drho = x.bk2.spacing;
Nrho = x.bk2.points;

rho_points = (0:(Nrho-1))*drho + rho0;

for i = 1:n
    x0(i,:) = interp1(rho_points,x.dat(i,:),rho);
end

