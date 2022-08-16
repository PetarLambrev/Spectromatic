function plotbygroup(Dat,Idx,GroupVars,varargin)
% PLOTBYGROUP - Plot spectra in Dat by groups defined in Idx by GroupVars
%
% Synthax
%   plotbygroup(Dat,Idx)
%   plotbygroup(Dat,Idx,GroupVars)
%   plotbygroup(Dat,Err,Idx,GroupVars)
%
% Description
%   plotbygroup(Dat,Idx,GroupVars) plots spectra in Dat in groups.
%   GroupVars is a list of variables in the table (categorical index) Idx.
%   
%   Shows the group variables in the title and the remaining variable in
%   the index table as legend text.
%
%   Uses the function plotfun for actual plotting. Modify plotfun to
%   control the plot appearance.
%
%   plotbygroup(Dat,Err,Idx,GroupVars) or
%   ploterrbygroup(Dat,Err,Idx,GroupVars) plots spectra and shaded errors.

% Input arguments
%   Dat (specdata) - spectra to plot 
%   Idx (table) - categorical index. Each row corresponds to one spectrum
%   GroupVars (string array) - grouping variables (must exist in Idx)
%
% See also: PLOTERRBYGROUP

    % Input validation
    if isa(Idx,'specdata') && istable(GroupVars)
        ploterrbygroup(Dat,Idx,GroupVars,varargin{:})
        return
    end
    if height(Idx) ~= numel(Dat)
        error('The number of rows in Idx must be equal to the number of spectra in Dat.')
    end
    if ~istable(Idx)
        error('Idx must be a table')
    end
   

    % Create title and legend strings from variable names
    
    if ~exist("GroupVars","var") || isempty(GroupVars)
        gi = true(1,width(Idx));
    else
        gi = contains(Idx.Properties.VariableNames,GroupVars);    
    end
    varstring = @(T) join(string(table2cell(T)),2);
    TitleText = varstring(Idx(:,gi));
    LegendText = varstring(Idx(:,~gi));

    % Plot
    G = findgroups(Idx(:,gi));
    if nargin > 3
        splitapply(@plotfun, Dat, TitleText, LegendText, varargin{:}, G)
    else
        splitapply(@plotfun, Dat, TitleText, LegendText, G)
    end
    
end
