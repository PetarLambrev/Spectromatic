function ploterrs(Dat,Err,Options)
% PLOTERRS - Plot spectra with errors in new figure
%
% Synthax
%   ploterrs(Dat,Err)
%   ploterrs(...,Name,Value)
%   h = ploterrs(...)
% 
% Description
%   ploterrs(Dat,Err) creates a new figure and plots Dat
%   ploterrs(...,Name,Value) allows for various formatting parameters
%   h = ploterrs(...) returns handle to the created figure
%
%  Note: ploterrs is primarily used in splitop and splitapply to create 
%     multiple figures
%
% Name-Value arguments
%     FigStyle - See FIGSTYLE
%     ColorScheme - See COLORSCHEME
%     Figure (struct) - Figure properties
%     TitleText - title text
%     LegendText - legend text
%     LegendFun - legend function
%     Axes (struct) - Axes properties
%     Legend (struct) - Legend
%     Line (struct) - Line properties
%     StartFun (function_handle) - pre-formatting function
%     EndFun (function_handle) - post-formatting function

arguments
    Dat specdata
    Err specdata
    Options.TitleText string = string.empty;
    Options.LegendText string = "ID";
    Options.LegendFun function_handle = function_handle.empty;
    Options.FigStyle = string.empty;
    Options.ColorScheme = string.empty;
    Options.Figure = struct.empty;
    Options.Axes = struct.empty;
    Options.Legend = struct.empty;
    Options.Line = struct.empty;
    Options.StartFun function_handle = function_handle.empty;
    Options.EndFun function_handle = function_handle.empty;
end

%% Create and pre-format figure
fh = figure;

% Run figureStart
if isempty(Options.StartFun) && exist("figureStart","file")==2
    figureStart;
end

fa = gca;
if ~isempty(Options.FigStyle)
    setfigstyle(Options.FigStyle) % note: External function
end
if ~isempty(Options.ColorScheme)
    setcolorscheme(ColorScheme) % note: External function
end

% Run Start function
if ~isempty(Options.StartFun)
    Options.StartFun(Dat)
end

%% Plot spectra
if isempty(Options.Line)
    LineArgs = {};
else
    LineArgs = namedargs2cell(Options.Line);
end

ploterror(Dat,Err,"LegendText",Options.LegendText,...
       "LegendFun",Options.LegendFun,LineArgs{:})
   
%% Post-formatting
if ~isempty(Options.TitleText)
    title(Options.TitleText(1))
end

% Set Figure parameters
if ~isempty(Options.Figure)
    set(fh,Options.Figure);
end

% Set Axes parameters
if ~isempty(Options.Axes)
    set(fa,Options.Axes);
end

% Set Legend parameters
if ~isempty(Options.Legend)
    fl = legend;
    set(fl,Options.Legend)
end

% Run End function
if ~isempty(Options.EndFun)
    Options.EndFun(Dat)
end

%% Return figure handle
if nargout
    h = fh;
end
end
