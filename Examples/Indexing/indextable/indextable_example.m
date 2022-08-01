% indextable_example
% 
% Example for creating categorical index with indextable. 
% This script performs the following types of operations:
%
% - Create a categorical index with indextable
% - Use splitop to average data across different groups
% - Use splitop to plot data per group

% Data contain CD and Abs for different mutants (WT, CM9, CM13),
% treated with different detergent (DM) concentrations.
% This information is NOT contained in the data but in an external file
% 'Table.xlsx'. The first column contains the file ID.

clear; close all

% Load Data (CD, Abs)
load data

% Load Excel table
tbl = readtable('Table.xlsx', 'ReadRowNames', true);
tbl.DM = tbl.DM * 100; % convert number to %

% Create index
C = indextable({CD.ID}, tbl);

% Average spectra using Mutant and DM as grouping variables
avg = CD.splitop(@mean, C, {'Mutant', 'DM'});
C = indextable({avg.ID}, tbl); % update index

% Sort by Mutant, DM
[C, ri] = sortrows(C, {'Mutant', 'DM'});
avg = avg(ri);

% Display index
disp(C)

% Plot Spectra per group

% First create custom legend text that indicates DM concentration
legfun = @(a) sprintf('DM %1.2f%%', a); % for example 'DM 0.07%'
legtxt = rowfun(legfun, C(:, 'DM'), 'OutputFormat', 'cell');

% Then use splitop to plot data grouped by Mutant
avg.set('ID', legtxt).splitop(@plotf, C, {'Mutant'}, [], 'PassGroupVars');