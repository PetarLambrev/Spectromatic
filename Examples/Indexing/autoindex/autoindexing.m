% Spectr-O-Matic Example - automatic indexing
%
% This is an example of using autoindex to group spectra 
% The data file contains a set of spectra for different samples
% There are four words in every spectrum ID representing grouping factors:
%   1 - aba4 or wt
%   2 - TM1 or TM2
%   3 - aDM or bDM
%   4 - fv or iso
%   The total number of spectra is 2^4 = 16.
%   In addition there are two baseline spectra - 'blank fv' and 'blank iso'
%
% The script first subtracts the relevant baseline from each spectrum,
% then aggregates TM1 and TM2 spectra together,
% and groups the iso and fv spectra in two figures.
%
% At the end, difference iso spectra aba4-wt are calculated.

clear; close all
load data

% Create autoindex x
x = dat.autoindex;

% x will have fields for every keyword (aba4, wt, TM1, TM2, ...)
% each field is a 1x18 logical array representing the keyword's occurrence

% show index map
disp(x)

% subtract blank fv and blank iso from all spectra 
dat(x.iso) = dat(x.iso) - dat(x.iso & x.blank);
dat(x.fv) = dat(x.fv) - dat(x.fv & x.blank);
dat(x.blank) = [];

% create a new index for newdat, which contains 16 spectra
x = dat.autoindex;

% assume that TM1 and TM2 are repetitions and average them
dat = (dat(x.TM1) + dat(x.TM2)) / 2;

% newdat has now only 8 (averaged) spectra, so recreate the index
x = dat.autoindex;

% plot the fv and iso spectra separately
figure; plot(dat(x.fv)); title('FV')
figure; plot(dat(x.iso)); title('ISO')

% calculate difference iso spectra aba4-wt
dif = dat(x.aba4 & x.iso) - dat(x.wt & x.iso);
dif = dif.set('YType', 'Difference CD');
figure; plot(dif); title('ISO aba4 - wt')