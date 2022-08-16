function ploterrfun(Dat,Err,TitleText,LegendText)
% PLOTERRFUN - Plot spectra with errors in new figure

figure;
ploterror(Dat,Err)
legend(LegendText)
title(TitleText(1))

end
