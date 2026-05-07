% ========================================================================
% SIGNAL PREPARATION: TRANSITION TIMES
%
% [x,t] = signal_transition_times(values,instants,ttrans,x0);
%
% Function that given a sequence of changes, and times, reconfigures the
% times in order to explicitly define when a change in a signal is starting
% and when its complete, so the interpolation later on keeps the sharp
% changes.
%
% 
% INPUTS
%
%   > values = array containing the values of the signal at each instant.
%   > instants = array with the times instants at which the changes occur
%                (it will contain also the end time, so it will have one
%                more element than the values array).
%   > ttrans = transition time for each change in the signal.
%   > x0 = initial value before the first change (0 by default)
%
% 
% OUTPUTS
%
%   > x = array containing the values reconfigured
%   > t = array containing the time instants reconfigured
%
% ========================================================================
% Hadriano Morales Escamilla
% 01/02/2021
% The University of Sheffield

function [x,t] = signal_transition_times(values,instants,ttrans,x0)

    % DEFAULT INITIAL VALUES IF NOT ENTERED
    if nargin < 4
        x0 = 0;
    end

    % VALUES AND INSTANS TO BE REDEFINE AS HORIZONTAL ARRAYS
    values = values(:)';
    instants = instants(:)';
    
    % TIME SEQUENCE
    t1 = instants(1:end-1); % when the changes start    
    t2 = t1 + ttrans; % when the changes finish
    t = [ t1 ; t2 ]; t = t(:)'; % combination of both
    
    % VALUES SEQUENCE
    x1 = [ x0 , values(1:end-1) ]; % values before the change is started
    x2 = values; % values after the change    
    x = [ x1 ; x2 ]; x = x(:)'; % combination of both
    
    % ADDITION OF LAST SAMPLE
    t = [ t instants(end) ]';
    x = [ x values(end) ]';
    
end