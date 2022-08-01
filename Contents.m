% Spectr-O-Matic Toolbox
% version 1.21
%
% This is a list of classes and functions in Spectr-O-Matic.
%
% Spectr-O-Matic classes
%
% specdata - one-way (X/Y) spectroscopy data
% specdata2D - two-way (e.g. time-resolved) data
% specparent - Common methods for 1D and 2D data
%
% Load data in Spectr-O-Matic
% 
% specdata/specdata - create specdata objects from MATLAB variables
% specdata/load - load data from text files
% specdata/unstack - unstack specdata array to specdata2D array
% specdata2D/specdata2D - create specdata2D objects from MATLAB variables
%
% Get and set properties
% 
% specdata/datatable, specdata/dt - MATLAB table with data and properties
% specdata/dim - number of data points in a spectrum
% specdata/get - get properties
% specdata/peaks - find local minima and maxima of Y
% specdata/proptable, specdata/pt - MATLAB table with properties
% specdata/set - set properties
% specdata/Yx - get Y value at given X
% specdata/xind - get indices of X values
% specdata/xymat - matrix of X/Y values
% specdata/xytable - table of X/Y values
% specdata2D/swapxt - exchange X and T axes (transpose Y)
% specdata2D/tind - get indices of T values
% specdata2D/Yt - get Y value at given T
% specdata2D/Yxt - get Y value at given X, T
%
% Mathematical operations
% 
% +, -, *, /, ^ - arithmetic operators
% specdata/mldivide - linear least squares fit
% specdata/abs - absolute values of Y
% specdata/baseline - subtract straight baseline
% specdata/diff - differentiate
% specdata/int - integrated area
% specdata/log, specdata/log10 - logarithm of Y
% specdata/min, specdata/max - minimum/maximum Y value
% specdata/mean - average spectra
% specdata/norm - normalize spectra
% specdata/smooth - smooth spectra
% specdata/sum - sum spectra
%
% Other operations on spectra
% 
% specdata/bin - bin data (boxcar) along X
% specdata/fit - fit model to spectra
% specdata/merge - join spectra into one
% specdata/setx - set new X axis and interpolate data
% specdata/setxlim - trim data to X range
% specdata/shiftx - shift spectra along X axis
% specdata/sort - sort spectra in an array based on a given property
% specdata2D/bint - bin data (boxcar) along T
% specdata2D/sett - set new T axis and interpolate data
% specdata2D/settlim - trim data to T range

%
% Search and indexing
% 
% specdata/autoindex - create a keyword index from spectra ID's
% specdata/catfind - search for keywords and return a categorical array
% specdata/catindex - create a catalog index by searching for keywords
% specdata/find - find spectra by keyword
% specdata/findindex, specdata/fi - find keywords and return a logical array
% indextable - match an index table to data
% specdata/splitop - split-apply operation on groups of spectra
% specdata/splitbinop - split-apply binary operation on groups of spectra
%
% Plotting
% 
% specdata/markpeaks - mark peak wavelengths on plot
% specdata/plot - plot spectra
% specdata2D/plot2D - surface plot of 2D spectra
% specdata2D/plotty - plot Y vs T
% specdata2D/plotxy  - plot Y vs X
%
% Save and export
% 
% specdata/save - save as a text file