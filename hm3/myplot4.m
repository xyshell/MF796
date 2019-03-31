function myplot4(X1, Y1, xlab)
%CREATEFIGURE(X1, Y1)
%  X1:  vector of x data
%  Y1:  vector of y data

%  Auto-generated by MATLAB on 10-Feb-2019 22:49:29

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create plot
plot(X1,Y1);

% Create ylabel
ylabel({'volatility'});

% Create xlabel
xlabel({xlab});

% Create title
title({'Implied Volatility Surface'});

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',13,'FontWeight','bold');
