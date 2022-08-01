function s = num2cellstr(n, formatSpec)
% num2cellstr - convert numeric array to cell array of strings
%
% Syntax
% s = num2cellstr(n)
% s = num2cellstr(n,precision)
% s = num2cellstr(n,formatSpec)
%
% Converts numeric array n to a cell array of strings s.
% Each element in n is represented as a string in the cell array s.
% precision - maximum number of significant digits
% formatSpec - format string
%
% see also - num2str

if ~isnumeric(n)
    error('Input must be a numeric array');
end

s = cell(size(n));

if ~exist('formatSpec', 'var')
    formatSpec = '%g';
end

for k = 1:numel(n)
    s{k} = num2str(n(k), formatSpec);
end
