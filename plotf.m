function h = plotf(x, GroupVars, FigStyle, ColorScheme)
% PLOTF - Plot spectra with a grouping variable title
%
% Synthax
%   plotf(x, GroupVars, FigStyle, ColorScheme)
%
% Example
%   data.splitop(@plotf, Idx, 'Group', 'PassGroupVars')
%
% Plots data split by the variable 'Group' of the categorical index Idx
% using the values of 'Group' as figure titles.

fighandle = figure;
plot(x)
if exist('GroupVars','var') && ~isempty(GroupVars)
    if iscell(GroupVars) && numel(GroupVars) > 1
        GroupVars = array2str(GroupVars);       
    end
    title(strjoin(GroupVars))
    legend boxoff
end

if exist('FigStyle','var')
    setfigstyle(FigStyle)
end

if exist('ColorScheme','var')
    setcolorscheme(ColorScheme)
end

if nargout
    h = fighandle;
end