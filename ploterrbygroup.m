function ploterrbygroup(Dat,Err,GroupVars,varargin)
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
%
% Name-Value arguments
%   See PLOTERRS for figure formatting Name-Value arguments

    % Input validation
    if numel(Err) ~= numel(Dat) 
        error('Dat and Err must have the same number of spectra.')
    end
    
    if nargin > 3
        argstruct = struct(varargin{:});
        if isfield(argstruct,'Index')        
            Idx = argstruct.Index;
            if ~istable(Idx)
                error('Idx must be a table')
            elseif height(Idx) ~= numel(Dat)
                error('The number of rows in Idx must be equal to the number of spectra in Dat.')
            end
            Dat = Dat.setmetadata(Idx);
            argstruct = rmfield(argstruct,'Index');
            varargin = namedargs2cell(argstruct);
        else
            Idx = Dat.proptable;
        end
    else
        Idx = Dat.proptable;
    end

    % Create title and legend strings from variable names
    gIndex = contains(Idx.Properties.VariableNames,GroupVars);    
    varstring = @(T) join(string(table2cell(T)),2);
    TitleText = varstring(Idx(:,gIndex));

    % Plot
    G = findgroups(Idx(:,gIndex));
    for gn = 1:max(G)
        gi = G==gn;
        ploterrs(Dat(gi),Err(gi),'TitleText',TitleText(gi),varargin{:})
    end
end