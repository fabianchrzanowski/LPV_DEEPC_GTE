% ========================================================================
% ARRAY LOOKUP
% 
% Generation of structure, or lookuptable object, to be used in simulink as
% lookup array.
%
% The case when instead of an array. it is onyl a value, is also handled
% here
%
% INPUTS
%   y: Array data, with two dimensions (m x p), where the array is 
%      described by the first dimension, and the second one is used for
%      scheduling
%   x: scheduling parameter values (as many elements as the second
%      dimension of y).
%   name: array name
%   lut_obj: option to indicate that a simulink lookuptable object is to be
%            created (false by default)
%
% OUTPUTS
%   LUT: structure containing the fields for the lookuptable implementation
%        in simulink, or the object (if he option is selected)
%
% ========================================================================
% Hadriano Morales Escamilla
% The University of Sheffield
% 22/07/2022

function LUT = lookup_array(y,x,name,lut_obj)
    
    % optional inputs
    if nargin < 4; lut_obj = []; end
    if nargin < 3; name = []; end
    
    % default option
    if isempty(lut_obj); lut_obj = false; end

    % default name
    if isempty(name); name = ''; end

    % number of elements
    m = size(y,1);

    % lookup table object creation (if selected), it assumes evenly spaced
    % scheduling space)
    if lut_obj
    
        % Preallocation
        LUT = Simulink.LookupTable;
        LUT.BreakpointsSpecification = 'Even spacing';

        % if more than one element in the array, the number of elements is
        % used as first breakpoint
        if m > 1

            % element
            LUT.Breakpoints(1).FieldName = 'element';
            LUT.Breakpoints(1).FirstPoint = 1;
            LUT.Breakpoints(1).Spacing = 1;
            LUT.Breakpoints(1).Min = 1;
            LUT.Breakpoints(1).Max = m;
    
            % scheduling variable
            LUT.Breakpoints(2).FieldName = 'sch';
            LUT.Breakpoints(2).FirstPoint = x(1);
            LUT.Breakpoints(2).Spacing = x(2)-x(1);
            LUT.Breakpoints(2).Min = x(1);
            LUT.Breakpoints(2).Max = x(end);

        % otherwise, there is only one breakpoint, with the scheduling
        % variable
        else

            % scheduling variable
            LUT.Breakpoints.FieldName = 'sch';
            LUT.Breakpoints.FirstPoint = x(1);
            LUT.Breakpoints.Spacing = x(2)-x(1);
            LUT.Breakpoints.Min = x(1);
            LUT.Breakpoints.Max = x(end);

        end

        % data
        LUT.Table.Value = y;
        LUT.Table.FieldName = name;

        % information for code generation
        LUT.StructTypeInfo.Name = ['LUT_array_' name];
        LUT.StructTypeInfo.HeaderFileName = 'header__lookups';

    % structure (if object not required)
    else

        % if more than one element in the array, the number of elements is
        % used as first breakpoint
%         if m > 1
            LUT.bk1.first = 1;
            LUT.bk1.spacing = 1;
            LUT.bk1.points = m;
            LUT.bk2.first = x(1);
            LUT.bk2.spacing = x(2)-x(1);
            LUT.bk2.points = length(x);
        % otherwise, there is only one breakpoint, with the scheduling
        % variable
%         else
%             LUT.bk1.first = x(1);
%             LUT.bk1.spacing = x(2)-x(1);
%             LUT.bk1.points = length(x);
%         end
        % data
        LUT.dat = y;
        % Name
        LUT.name = ['LUT_array_' name];

    end


end

