% ========================================================================
% MATRIX LOOKUP
% 
% Generation of structure, or lookuptable object, to be used in simulink as
% lookup table.
%
% INPUTS
%   M: Matrix data, with three deimensions (m x n x p), where the matrix is
%       decribed by the two first dimensions, and the third one is used for
%       scheduling
%   x: scheduling parameter values (as many elements as the thrid dimension
%       M).
%   name: matrix name
%   lut_obj: option to indicate that a simulink lookuptable object is to be
%            created (flase by default)
%
% OUTPUTS
%   LUT: structure containing the fields for the lookuptable implementation
%           in simulink, or the object (if he option is selected)
%
% ========================================================================
% Hadriano Morales Escamilla
% The University of Sheffield
% 22/07/2022

function LUT = lookup_matrix(M,x,name,lut_obj)
    
    % optional inputs
    if nargin < 4; lut_obj = []; end
    if nargin < 3; name = []; end
    
    % default option
    if isempty(lut_obj); lut_obj = false; end

    % default name
    if isempty(name); name = ''; end


    % rows and columns
    [m,n,~] = size(M);

    % lookup table object creation (if selected), it assumes evenly spaced
    % scheduling space)
    if lut_obj
    
        % Preallocation
        LUTObj = Simulink.LookupTable;
        LUTObj.BreakpointsSpecification = 'Even spacing';

        % rows
        LUTObj.Breakpoints(1).FieldName = 'row';
        LUTObj.Breakpoints(1).FirstPoint = 1;
        LUTObj.Breakpoints(1).Spacing = 1;
        LUTObj.Breakpoints(1).Min = 1;
        LUTObj.Breakpoints(1).Max = m;
        
        % columns
        LUTObj.Breakpoints(2).FieldName = 'column';
        LUTObj.Breakpoints(2).FirstPoint = 1;
        LUTObj.Breakpoints(2).Spacing = 1;
        LUTObj.Breakpoints(2).Min = 1;
        LUTObj.Breakpoints(2).Max = n;

        % scheduling variable
        LUTObj.Breakpoints(3).FieldName = 'sch';
        LUTObj.Breakpoints(3).FirstPoint = x(1);
        LUTObj.Breakpoints(3).Spacing = x(2)-x(1);
        LUTObj.Breakpoints(3).Min = x(1);
        LUTObj.Breakpoints(3).Max = x(end);

        % data
        LUTObj.Table.Value = M;
        LUTObj.Table.FieldName = name;

        % information for code generation
        LUTObj.StructTypeInfo.Name = ['LUT_matrix_' name];
        LUTObj.StructTypeInfo.HeaderFileName = 'header__lookups';

    % structure (if object not required)
    else

        LUT.bk1.first = 1;
        LUT.bk1.spacing = 1;
        LUT.bk1.points = m;
        LUT.bk2.first = 1;
        LUT.bk2.spacing = 1;
        LUT.bk2.points = n;
        LUT.bk3.first = x(1);
        LUT.bk3.spacing = x(2)-x(1);
        LUT.bk3.points = length(x);
        LUT.dat = M;
        LUT.name = ['LUT_matrix_' name];

    end


end

