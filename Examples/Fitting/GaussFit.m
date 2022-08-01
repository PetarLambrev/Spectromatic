% SpecData example: gaussian fit
% 
% This script uses the SpecData class 
% to fit spectral data using gaussian model
%
% Petar Lambrev, 2011

clear; close all;                   % clear memory
 
%% Load data from ASCII files
Abs = specdata.load('abswt.txt');
Abs = Abs.setxlim([600 800]);

%% Fit using 4-gaussian model
fitres = Abs.fit('gauss4');

% display result
disp(fitres)

% compute individual components
cmp = gaussians(fitres, Abs.X);

%% Plot result
figure;

% plot fit result
plot(fitres,Abs.X,Abs.Y);
hold on

% overlay individual components
cmp.plot('Color','r');
legend off