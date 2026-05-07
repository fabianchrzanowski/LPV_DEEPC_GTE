% AGTF30 engine setup file
% ========================================================================
% 
%
% Pablo Baldivieso Monasterios
% The University of Sheffield
% 18/12/2023

MWS.engName = 'AGTF30';
MWS = setup_Controller(MWS);
MWS = setup_AllEng(MWS);

% load bus
evalin('base','load(''Eng_Bus.mat'')');