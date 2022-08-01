% Spectr-O-Matic Example - Spectral decomposition
%
% This example decomposes a mixture spectrum into reference components
% by linear least squares fit (left matrix division).

clear; close all

%% Load reference and mixture spectra
% A - 2 gaussian spectra
% B - a 3:2 mixture spectrum B = 3*A(1) + 2*A(2)
load mixturedata

figure; plot(A) % reference spectra
figure; plot(B) % spectrum of the mixture

%% Decompose mixture spectrum
c = A \ B;
disp(c)

% Inspect fit
F = sum(A*c); F.ID = 'Fit';
R = F - B;    R.ID = 'Residuals';
figure; plot(B, F, R)
