% Plotting smooth curves example

clear; close all;

% load data
Abs = specdata.load('abswt.txt');

% bin every 30 points
Abs = Abs.bin(50);

% plot binned data
Abs.plot('Marker','o');

% plot data using smoothing
Abs.plot('SmoothLine','Color','r');

legend('binned', 'smooth');