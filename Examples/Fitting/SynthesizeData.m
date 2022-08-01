% Spectr-O-Matic Example - Synthesize data for linear decomposition
%
% This example creates a synthetic dataset for use in the 
% linear decomposition example

clear; close all

% Create 2 reference spectra as gaussians
x = (0:100)'; % X axis values
gauss = @ (a, b, c, x) a.*exp(-(x-b).^2/(2*c^2)); % gaussian function
y1 = gauss(0.5, 30, 20, x); % 1st gaussian 
y2 = gauss(1.0, 60, 10, x); % 2nd gaussian

A = specdata(x, [y1,y2], {'A1','A2'}); % reference spectra
figure; plot(A)

% Create a "mixture" spectrum
c1 = 3; c2  = 2;                       % mixture coefficients
y = c1.*y1 + c2.*y2;                   % create mixture
y = y + 0.2*rand(numel(x),1) - 0.1;    % add random noise

B = specdata(x, y, 'Mixture');         % create spectrum
figure; plot(B)

% Save Data
save mixturedata A B
