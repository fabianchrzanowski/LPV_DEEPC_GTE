%% LPV simulation
% ------------------------------------------------------------------------
% Simlate LPV system with LPV_sys a structure containing the matrices
% defining the LPV system and their domain of definition. (u,d) are vectors 
% of inputs and disturbances, and rho is the operating point. 
% 
% the ouput of the procedure is the state and output sequences.
%
%
% Pablo R Baldivieso Monasterios
% 14/08/2022
% The University of Sheffield
% ------------------------------------------------------------------------
function [y,x,e] = LPVsim(LPV_sys,y_uf,u,d,rho)

T = size(u,2);
p = size(LPV_sys.C,1);
n = size(LPV_sys.A,1);
m = size(LPV_sys.B,2);

y = zeros(p,T);
e = zeros(p,T);
x = zeros(n,T);

rho0 = rho(:,1);
ud  = [u;d];

for tt = 1:T

    if tt == 1 || sum(rho(:,tt) == rho0) < 3
        A = LPV2mat(LPV_sys.A,LPV_sys.rho,rho(:,tt));
        B = LPV2mat(LPV_sys.B,LPV_sys.rho,rho(:,tt));
        C = LPV2mat(LPV_sys.C,LPV_sys.rho,rho(:,tt));
        D = LPV2mat(LPV_sys.D,LPV_sys.rho,rho(:,tt));
        rho0 = rho(:,tt);
    end


    if tt == 1
        % initial state
        
        H = [      D       zeros(p,6*m);...
                 C*B       D       zeros(p,5*m);...
               C*A*B     C*B       D       zeros(p,4*m);...
             C*A^2*B   C*A*B     C*B       D     zeros(p,3*m);...
             C*A^3*B C*A^2*B   C*A*B     C*B     D   zeros(p,2*m);...
             C*A^4*B C*A^3*B C*A^2*B   C*A*B   C*B   D zeros(p,m);...
             C*A^5*B C*A^4*B C*A^3*B C*A^2*B C*A*B C*B D];
    
        O  = [C;C*A;C*A^2;C*A^3;C*A^4;C*A^5;C*A^6];
        x(:,1) = linsolve(O,(vec(y_uf(:,1:7)) - H*vec(ud(:,1:7))));
    end
    
    x(:,tt+1) = A*x(:,tt) + B*ud(:,tt);
    y(:,tt)   = C*x(:,tt) + D*ud(:,tt);
    e(:,tt)   = (y_uf(:,tt) - y(:,tt))./y_uf(:,tt);
end