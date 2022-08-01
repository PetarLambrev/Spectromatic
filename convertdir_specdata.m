% Batch convert to specdata

clear

% Find files *.mat in all subfolders
[~, files] = system('dir /s/b *.mat');
files = textscan(files, '%s', 'Delimiter', '\n');
files = files{1};

% convert all
for k = 1:numel(files)
    convert_specdata(files{k});
end
