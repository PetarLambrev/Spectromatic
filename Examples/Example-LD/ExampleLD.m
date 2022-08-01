% Spectr-O-Matic Example
% Analysis of linear dichroism data
% 
% This script uses the specdata class 
% to manipulate linear dichroism data.
%
% Petar Lambrev, 2017

clear; close all;                   % clear memory
 
%% Load and index data

% Load data from text files into a specdata object
data = specdata.load('Data\LD*.txt');
Abs = specdata.load('Data\Abs.txt');

%% Perform calculations
% Get only first column of data
data = data(:,1);

% Subtract the blank reference spectra
data = data.find('Chloroplast') - data.find('Blank');

% Create a new keyword index for the LD data
x = data.autoindex;

% Subtract vertical and horizontal spectra
LD(1) = data(x.EH) - data(x.EV); % edge-aligned sample
LD(2) = data(x.FH) - data(x.FV); % face-aligned sample

% Create new ID
LD = LD.set('ID', {'edge-aligned', 'face-aligned'});

% Bin every 5 datapoints together
LD = LD.bin(5);          

% Limit the wavelength range to 400-750
LD = LD.setxlim([400 750]);    

% Smooth spectra
LD = LD.smooth;                

% Shift Abs spectrum to cross zero at 800 nm
Abs = Abs - Abs.Yx(800);

% Set X axis of Abs to match LD
Abs = Abs.setx(LD(1).X);

% Calculate Anisotropy as LD/Abs
r = LD / Abs;
r = r.set('YType','LD/A');

%% Plot results

% Plot LD
figure;                         % create a figure
plot(LD);                       % plot the spectrum
legend('location','northwest'); % change legend position
markpeaks(LD,0.01);             % mark peak positions on the graph

% Plot anisotropy
figure;
r.plot;