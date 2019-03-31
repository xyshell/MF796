function myplot1(X1, Y1, X2, Y2)
%CREATEFIGURE(X1, Y1, X2, Y2)
%  X1:  vector of x data
%  Y1:  vector of y data
%  X2:  vector of x data
%  Y2:  vector of y data

%  Auto-generated by MATLAB on 27-Feb-2019 21:54:37

% Create figure
figure1 = figure('NumberTitle','off','Name','LiveEditorFigure');

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create plot
plot(X1,Y1,'DisplayName','1M');

% Create plot
plot(X2,Y2,'DisplayName','3M');

% Create ylabel
ylabel({'vol'},'FontWeight','bold');

% Create xlabel
xlabel({'K'});

% Create title
title({'Implied Vol Surface'});

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',13,'FontWeight','bold');
% Create legend
legend(axes1,'show');

