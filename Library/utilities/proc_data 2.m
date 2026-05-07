function data = proc_data(z,names,t,data,f)
% ========================================================================
% Process data acording to min-max such that each signal x\in[-1,1]
% 
% ========================================================================
%
%
% Pablo Baldivieso Monasterios
% The University of Sheffield
% 01/07/2024

nz = numel(names);

for iz = 1:nz
    z(iz) = z(iz).resample(t);
    z(iz).Name = names{iz};

    data.offset(iz,f)     = z(iz).Data(1);
    data.scale.max(iz,f)  = 1;%max(abs(z(iz).Data - data.offset(iz,f)));
    
    data.raw(:,iz,f)      = z(iz).Data;
    data.proc(:,iz,f)     = (z(iz).Data - data.offset(iz,f))/data.scale.max(iz,f);
%   data.scale.min(iz,f)  = min((z(iz).Data));
%   data.proc(:,iz,f)     = (2*z(iz).Data - data.scale.max(iz,f) - data.scale.min(iz,f))/...
%                              (data.scale.max(iz,f) - data.scale.min(iz,f));
end

end

