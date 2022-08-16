function h = plotfun(Dat,TitleText,LegendText,FigStyle,ColorScheme)
% PLOTFUN - Plot spectra with errors in new figure

fh = figure;
plot(Dat)
legend(LegendText)
title(TitleText(1))

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
