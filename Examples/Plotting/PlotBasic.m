% plot spectra example
% 
% This example creates an array of two SpecData objects
% containing sine and cosine function
% and plots them on a graph

x = 0:2:360; ysin = sind(x); ycos = cosd(x);
l = struct('XType','Angle','XUnit','deg','YType','value');

trig(1) = specdata(x, ysin, 'sin', l);
trig(2) = specdata(x, ycos, 'cos', l);

figure; trig.plot;