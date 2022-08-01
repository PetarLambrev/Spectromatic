%% Spectr-O-Matic example: catindex
% Example script for Spectr-O-matic toolbox. 
% 
% This script performs the following types of operations:
% 
% * Load and combine multiple datasets
% * Create a categorical index with *catindex*
% * Use *splitbinop* to subtract control sample for different groups
% * Use *splitop* to average data across different groups
% * Use splitop to plot data per group
% * Create a pivot table for a calculated parameter
% 
% Data description:
% 
% * Spectra from two experiments are stored as text files in folders _Exp1_, 
% _Exp2_
% * Text file names contain information about the sample: mutant (wt, cm9 or 
% cm13)  and treatment (nothing, 0.07% bdm, 0.1% bdm)
% * Every text file contains three spectra - CD, HT, Abs.
% 
% Brief algorithm:
% 
% # Load data
% # Subtract reference (baseline) spectra across experiments
% # Perform simple calculations
% # Create index - categorize data based on mutant and treatment
% # Perform calculations on groups of data
% # Plot spectra per group
% # Create a parameter table

clear; close all

%% 1) Load Data

% Load text files from Exp1 and Exp2
exp1 = specdata.load('Exp1\*.txt', 'FileType', 'JWS'); 
exp2 = specdata.load('Exp2\*.txt', 'FileType', 'JWS');
% Create combined dataset
dat = [exp1; exp2];

% Split CD and absorption spectra
CD = dat(:,1); Abs = dat(:,3);

%% 2) Subtract baseline across experiments
% # Create an index of reference spectra (containing the keyword _blank_)
% # Create a categorical index table |*T|* with a single variable _ExpID_
% # Subtract reference spectra from the rest of the spectra, using ExpID as 
% a grouping variable

bl  = CD.findindex('blank');  % find reference spectra
T   = CD.pt('ExpID');          % create a table T with one column containing ExpID
CD  = CD.splitbinop(@minus, ~bl, bl, T, 'ExpID')'; 
Abs = Abs.splitbinop(@minus, ~bl, bl, T, 'ExpID')';

%% 3) Create index
% Define variables

v.mutant = {'wt', 'cm9', 'cm13'};
v.dm = {'?0%', '0.07%', '0.1%'}; % default value is 0% (marked with ?)

% define fancy variable labels
w.mutant = {'WT', 'CM9', 'CM13'};

% Create index
C = CD.catindex(v, w);

% Add experiment ID (ExpID) as a variable in the index
C.ExpID = categorical({CD.ExpID}');

%% 4) Simple calculations
% Shift absorption to 0 at 750 nm

Abs = Abs - Abs.Yx(750);

% Normalize CD to absorption maximum between 600 and 700 nm
CD = CD / Abs.max([600 700]);

% Average spectra using mutant and dm as grouping variables
Avg = CD.splitop(@mean, C, {'mutant', 'dm'});
C = Avg.catindex(v, w); % update index

%% 5) Sort data
% Sort data by mutant, dm Note that categorical variables are sorted in the 
% order they are defined in the struct v

[C, ri] = sortrows(C, {'mutant', 'dm'});
Avg = Avg(ri);

%% 6) Plot Spectra per group
% # Plot averaged difference spectra grouped by mutant 
% # Write the mutant name in the figure title and dm concentration in the legend

dm = cellstr(C.dm); % dm concentration for each spectrum;
Avg.set('ID', dm).splitop(@plotf, C, {'mutant'}, [], 'PassGroupVars');

%% 7) Create a parameter table 
% # Calculate a parameter from the spectra and add it to the index table
% # Make a parameter table separating mutants into columns
% # Plot the table data

% Add calculated PeakRatio parameter
C.PeakRatio = Avg.Yx(665) ./ Avg.Yx(681);
% Make a wide table of PeakRatio with mutant in columns
PR = unstack(C,'PeakRatio','mutant');
disp(PR)