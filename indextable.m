function [sorted, vars] = indextable(keys, vartable, keyvar)
% Sort variable table using a list of keys
% 
% sorted = indextable(keys, vartable)
% sorted = indextable(keys, vartable, keyvar)
% [sorted, vars] = indextable(...)
%
% sorted = indextable(keys, vartable)
%
% Sort the rows in "vartable" matching row names to the cell array "keys".
% If vartable is a string, then the table is first read from that file.
%
% sorted = indextable(keys, vartable, keyvar)
%
% Sort the rows in "vartable" so that the column with name "keyvar" matches
% the cell array "keys"
%
% [sorted, vars] = indextable(...)
% 
% Return also a struct "vars" containing the unique values of all variables


%% Load vartable from file
if ischar(vartable)
    if exist('keyvar','var')
        vartable = readtable(vartable,'ReadRowNames',false);
    else
        vartable = readtable(vartable,'ReadRowNames',true);
    end
end
    
if ~exist('keyvar','var')
    keyvar = '';
end

% Convert cell arrays to categorical arrays
vars = vartable.Properties.VariableNames;
for k = 1:numel(vars)
    var = vars{k};
    v = vartable.(k);
    if iscell(v) && ~strcmp(var,keyvar) && isempty(regexpi(var,'(remark|comment)'))
        vartable.(k) = categorical(v);
    end
end

% Get right keys from vartable
if exist('keyvar','var') && ~isempty(keyvar)
    rightkeys = vartable.(keyvar);
else
    rightkeys = vartable.Properties.RowNames;
end

 %% find matching keys
[ia,ib] = ismember(lower(keys),lower(rightkeys));

if any(ia==0)
    missinglist = keys(ia==0);
    missingstring = '';
    for k = 1:numel(missinglist), missingstring = [missingstring, ' ', missinglist{k}]; end
    warning('Some keys were not found in the source index table: \n%s', missingstring)
    ib(ib==0) = [];
end

%% create sorted table
sorted = vartable(ib,:);
vars = extractvars(sorted);



