function h = plotfun(Dat,TitleText,LegendText,FigStyle,ColorScheme)
% PLOTFUN - Plot spectra with errors in new figure

fh = figure;
plot(Dat)
if ~isempty(LegendText)
    legend(LegendText)
end
if ~isempty(TitleText)
    title(TitleText(1))
end

if exist('FigStyle','var')
    setfigstyle(FigStyle)
end

if exist('ColorScheme','var')
    setcolorscheme(ColorScheme)
end

if nargout
    h = fh;
end

end
