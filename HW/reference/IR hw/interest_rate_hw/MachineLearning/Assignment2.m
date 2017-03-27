% script:  box_and_ema_compare.m
% descrip: hw3 solutions that study Box and Ema Comparison
% author:  Ray He

% Critical three lines
% ----------------------------------------------------------------------- %
close all           % close all open figure
clear variables
clc
% ----------------------------------------------------------------------- %

upbound = 1;

plotPerRow = 3;

data = [[0, 0]; [0, 1]; [1, 1]; [0, 2]; [1, 2]; [2, 2]; [5, 10]; ...
    [50, 100]; [500, 1000]];

rows = 3;

x2 = linspace(0,upbound);
y2 = 1 + x2 - x2;

std = 0.1;
mean = 0.5;

y2 = 1 / sqrt(2 * pi * std) * exp(-1 * (x2 - mean).^2 / (2 * std^2));

for i = 1 : length(data)
    N = data(i, 2);
    R = data(i, 1);

    y1 = x2.^R .* (1 - x2).^(N - R); 

    subplot(rows,plotPerRow,i);
 
    str = sprintf('%d heads out of %d tosses', R, N);
    
    [AX,H1,H2] = plotyy(x2,y1,x2,y2);

    title(str);    
    
    set(H2,'color','red')

    set(AX(2), 'ycolor', 'r');
    set(AX(2),'YLim',[0 2]);
    set(AX(2),'YTick',0:.5:2);
    
    ylabel(AX(2),'Prior  P(H|I)') % right y-axis
    ylabel(AX(1),'Posterior P(H|data, I)') % left y-axis

end