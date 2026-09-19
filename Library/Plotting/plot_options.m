%% Plotting options
% ------------------------------------------------------------------------
% Color vector for plotting (6 different options) in RGB coordinates
% default values for figures, i.e., latex
% interpreter, units, position, and font sizes
%
% Pablo R Baldivieso Monasterios
% 11/04/2022
% The University of Sheffield
% ------------------------------------------------------------------------

    color_plot = [[ 28 123 243];... red
                  [237 141  61];... blue green   
                  [128  83 180];... purple
                  [ 59 167  37];... light green
                    0 204   0;... green
                  204 102   0;... orange
                    0 204 204;... dark cyan
                    0 102 204;... navy
                    0   0 204;... blue
                  204 204 204;... yellow
                  204   0 204;... pink
                  204   0 102;... magenta
                   96  96  96 ... gray
                  ]/255;

    marker_plot = {'none','+','*','x','s','d','p'};
    % restore default
    reset(groot)

    
    % figure font size
    set(0,'defaultAxesFontSize', 12);

    % default color
    set(0, 'defaultFigurecolor',[1 1 1]);

    % default linewidth
    set(groot, 'DefaultLineLineWidth' , 2);
    set(groot, 'DefaultStairLineWidth', 2);

    % set interpreter to latex for axis, titles, legends, ticks
    set(groot, 'defaultTextinterpreter','latex');  
    set(groot, 'defaultAxesTickLabelInterpreter','latex');  
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot, 'defaultColorbarTickLabelInterpreter','latex');
    set(groot, 'defaultAxesXGrid','on')
    set(groot, 'defaultAxesYGrid','on')

    % Keep it simple — no global figure sizing or positioning to prevent crashes.