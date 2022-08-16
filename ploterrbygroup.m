function ploterrbygroup(Dat,Err,Idx,GroupVars,varargin)
% PLOTERRBYGROUP - Plot spectra in Dat by groups defined in Idx by GroupVars
%
% Synthax
%   ploterrbygroup(Dat,Err,Idx,GroupVars)
%
% Description
%   ploterrbygroup(Dat,Err,Idx,GroupVars) plots that spectra in Dat 
%   with shaded errors in Err separated in groups defined by GroupVars.
%   GroupVars is a list of variables in the table (categorical index) Idx.
%   
%   Shows the group variables in the title and the remaining variable in
%   the index table as legend text.
%
%   Uses the function ploterrfun for actual plotting. Modify ploterrfun to
%   control the plot appearance.
%
% Input arguments
%   Dat (specdata) - spectra to plot 
%   Err (specdata) - errors (variance)
%   Idx (table) - categorical index. Each row corresponds to one spectrum
%   GroupVars (string array) - grouping variables (must exist in Idx)

    % Input validation
    if numel(Err) ~= numel(Dat) 
        error('Dat and Err must have the same number of spectra.')
    end
    if height(Idx) ~= numel(Dat)
        error('The number of rows in Idx must be equal to the number of spectra in Dat.')
    end
    if ~all(contains(GroupVars,Idx.Properties.VariableNames))
        error('Variable not found. GroupVars must contain variable names found in Idx.')
    end
    if ~istable(Idx)
        error('Idx must be a table')
    end    

    % Create title and legend strings from variable names
    gi = contains(Idx.Properties.VariableNames,GroupVars);    
    varstring = @(T) join(string(table2cell(T)),2);
    TitleText = varstring(Idx(:,gi));
    LegendText = varstring(Idx(:,~gi));

    % Plot
    G = findgroups(Idx(:,gi));
    if nargin > 3
        splitapply(@ploterrfun, Dat, Err, TitleText, LegendText, varargin{:}, G)
    else
        splitapply(@ploterrfun, Dat, Err, TitleText, LegendText, G)
    end

end
