% Spectr-O-Matic Example - multiple baselines
%
% Example for subtracting multiple baselines from groups of data
% using group_baseline_corr
%
% In this example dat contain 9 spectra:
% 6 sample measurements, 3 baseline measurements
%
% 3 measurements per sample and 3 baselines - fh, fv, iso.
%
% The script subtracts the respective baseline from each spectrum,
% then plots the spectra.

clear; close all
load data

% Subtract baselines
dat2 = group_baseline_corr(dat, {'fh','fv','iso'}, 'blank');

% Display ID and History of operations
disp(dat2.pt({'ID', 'History'}))

% Plot corrected spectra
figure('name','corrected');
plot(dat2)
