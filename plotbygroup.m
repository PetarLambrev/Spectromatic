function plotbygroup(Dat,GroupVars,varargin)
% PLOTBYGROUP - Plot spectra in Dat by groups defined in Idx by GroupVars
%
% Synthax
%   plotbygroup(Dat,GroupVars)
%   plotbygroup(Dat,GroupVars,'Index',Idx)
%   plotbygroup(...,Name,Value)
%
% Description
%   plotbygroup(Dat,GroupVars) plots spectra in Dat in groups specified by
%   the list of property/metadata variables GroupVars.
%
%   plotbygroup(Dat,GroupVars,'Index',Idx) uses an external index table Idx.
%   
%   plotbygroup(...,Name,Value) allows various formatting parameters.
%
%   Uses the function plots for actual plotting. See PLOTS for details.
%
%   plotbygroup(Dat,Err,...) or
%   ploterrbygroup(Dat,Err,...) plots spectra and shaded errors.
%
% Input arguments
%   Dat (specdata) - spectra to plot 
%   GroupVars (string array) - grouping variables (must exist in Idx)
%
% Name-Value arguments
%   Index - External categorical index to use
%   See PLOTS for figure formatting Name-Value arguments
%
% See also: PLOTS, PLOTERRBYGROUP

    % Input validation
    if isa(GroupVars,'specdata')
        ploterrbygroup(Dat,GroupVars,varargin{:})
        return
    end
    
    if nargin > 2
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
        argstruct = struct;
        Idx = Dat.proptable;
    end

    % Create title and legend strings from variable names
    if isempty(GroupVars)
        gIndex = true(1,width(Idx));
    else
        gIndex = matches(Idx.Properties.VariableNames,GroupVars);    
        if ~any(gIndex)
            error('The grouping variables were not found in index/metadata');
        end
    end
    varstring = @(T) join(string(table2cell(T)),2);
    varnames = string(Idx.Properties.VariableNames);
    if isfield(argstruct,'TitleText') && matches(string(argstruct.TitleText),varnames,'IgnoreCase',true)
        TitleText = varstring(Idx(:,argstruct.TitleText));
        argstruct = rmfield(argstruct,'TitleText');
        varargin = namedargs2cell(argstruct);
    else        
        TitleText = varstring(Idx(:,gIndex));
    end

    % Plot
    G = findgroups(Idx(:,gIndex));
    for gn = 1:max(G)
        gi = G==gn;
        plots(Dat(gi), "TitleText", TitleText(gi), varargin{:})
    end
end
