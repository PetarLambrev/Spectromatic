%% Create specdata objects
% Create specdata objects from variables and text files

% Create some X and Y arrays
x = 0:0.1:pi;
y = sin(x);

% Create a spectrum
S = specdata(x,y,'sinx');

% Create a second spectrum
y = cos(x);
S(2) = specdata(x,y,'cosx');

% Plot the spectra
figure; plot(S)

