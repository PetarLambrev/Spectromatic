% Spectr-O-Matic example library
%
% This is a library of common Spectr-O-matic commands.
%
% How to use:
% 1. From the Editor tab, click Go To and select a section
% 2. Copy a line and paste into your script
% 3. Modify parameters, etc. as needed

return

%% Load Data

% Load *.txt files
CD = SpecData.load('*.txt');

% Load Jasco *.txt files
CD = SpecData.load('*.txt', 'FileType', 'JWS');

% Load Jasco *.txt files from directory 'Data'
CD = SpecData.load('Data\*.txt', 'FileType', 'JWS');

% Save object "CD" as file "data.mat"
save data CD

% Load data from file "data.mat"
load data

% Load data from file "data.mat" to object "dat1"
dat1 = load('data');

%% Modify properties

% Set properties
CD = CD.set('XType', 'Wavelength');
CD = CD.set('YType', 'CD');
CD = CD.set('YUnit', '\DeltaOD');
CD = CD.set('YType','CD','YUnit','mdeg');

%% Parse ID text

% Get spectra ID's from dat into variable IDs
IDs = {CD.ID}';

% Give new ID's to spectra
CD = CD.set('ID', IDs);

% Replace text in spectra IDs
CD = CD.set('ID', strrep({CD.ID}, 'oldtext', 'newtext'));

% Remove file extension from IDs
CD = CD.remove_ext;

% Remove -1, -2, -3, ... from the end of the ID
CD = CD.set('ID',regexprep({CD.ID},'-\d$',''));

% Read numeric data from ID
Num = str2double(regexp(Txt, '(?<=Temperature)(\d+)','match','once'));
% ID = {'Temperature15', 'Temperature21.75'} -> Num = [15, 21]

%% Subtract baseline

% Single baseline file
bl = CD.fi('baseline');
CD = CD(~bl) - CD(bl);

% Single baseline file, multiple columns
bl = CD(:,1).fi('baseline');
CD = CD(~bl,:) - CD(bl,:);

% Different baselines per group, multiple columns
% (Groups are iso, fh, ev, buffer. Single baseline for each group.)
dat = group_baseline_corr(dat, {'iso', 'fh', 'ev', 'buffer'}, 'baseline');

% Different baselines per group, one column, using index x
bl = x.Sample=='baseline';
CD = CD.splitbinop(@minus, ~bl, bl, x, {'Species','Treatment'})';

%% Trim and smooth
% Set X range
CD = CD.setXlim([350 750]);

% Smooth
Abs = Abs.smooth(5);
CD = CD.smooth(4, 'sgolay');

%% Shift and Normalize
% Shift to zero at a certain wavelength
CD = CD - CD.Yx(750);

% Normalize by maximum
dat = dat.norm;

% Normalize by maximum found in X range
dat = dat / dat.max([600 800]);

%% Calculate parameters
d = dat.Yx(515) - dat.X(535);
a = dat.max([600 800]) - dat.min([600 800]);
FvFm = 1 - (F.Yx(0.02) / F.max([500 2000]));
Farea = F.norm.int([0 1000]);

%% Create Index
% Create index by keywords in ID
x = CD.autoindex;

% Define variables
vars.Species = {'pea', 'spinach', 'Arabidopsis'};
vars.Mutant = {'wt', 'npq1', 'npq2'};

% Define fancy variable values
newv.Species = {'P. sativum', 'S. oleracea', 'A. thaliana'};

% Create index from variable list
x = CD.catindex(vars);

% Create index with fancy variable values
x = CD.catindex(vars, newv);

% Create index using variables in Excel table
x = indextable({CD.ID}, 'Workbook.xlsx')

% Create index using preloaded Excel table
wbk = readtable('Workbook.xlsx','ReadRowNames',true);
x = indextable({CD.ID}, wbk)

%% Sort data by category
[catx, rowi] = sortrows(catx, {'Var1', 'Var2'});
dat = dat(rowi);


%% Calculations using index

% Subtract reference spectra per group
dif = dat.splitbinop(@minus, catx.Sample~='control', catx.Sample=='control', catx, {'GroupVar1', 'GroupVar2'})';

%% Plot Data

% Simple plot
figure;
dat.plot;

%% Plot Data using index

% Plot by groups
CD.splitop(@plotf, x);

% Plot by groups with selected grouping variables
CD.splitop(@plotf, x, {'Species','Sample'});
CD.splitop(@plotf, x, {'Species','Sample'}, [], 'PassGroupVars');

% Create custom legend text from index variables
% using arrays2tr
legends = array2str({x.Medium, {CD.ID}});
legends = array2str({x.Medium, x.FileNo}, ' #');

% Create custom legend text with sprintf
legfun = @(a1,a2) sprintf('%s #%d', char(a1), a2);
newID = rowfun(legfun, x(:, {'Medium','FileNo'}), 'OutputFormat', 'cell');

% Plot by groups using custom legend text
CD.set('ID', legends).splitop(@plotf, x);
CD.set('ID', legends).splitop(@plotf, x, {'Species','Sample'}, [], 'PassGroupVars');

% Plot by groups using custom function
g = findgroups(x.Sample);                       % group data by Sample
splitapply(@plt, CD, x.Species, x.Sample, g)    % plot data groups using plt

function plt(d, t, l)
    % d - data, t - title, l - legend
    figure; plot(d); 
    title(sprintf('Sample: %s', t(1)))
    legend(cellstr(l))
end

