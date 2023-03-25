Spectr-O-Matic Toolbox for MATLAB(R)

by Petar H. Lambrev (2011-2022)

Version 2.4.1

Spectr-O-Matic is a MATLAB toolbox for analysis of spectroscopy data.  

Requires MATLAB R2020b or later.

For introduction and tutorials, see
https://plambrev.wixsite.com/spectromatic

Download the latest version at
http://www.mathworks.com/matlabcentral/fileexchange/32828-spectr-o-matic

Release notes

v. 2.4.2
- altmeta - add alternative human-friendly versions of metadata

v. 2.4.1

Main features

Spectra now can contain custom metadata. The custom metadata can be used 
for grouped operations, like averaging, plotting, etc.

The specparent class now enforces the data types of its properties, e.g. ID is string.

New methods

- metaindex creates a categorical index (catindex) and stores it as metadata
- setmetadata (setmd) assigns the contents of a table as metadata
- addmetadata (addmd) adds a custom metadata field to spectra
- deletemetadata (deletemd) deletes metadata from spectra
- metatable (mt) returns the custom metadata for all spectra as a table
- table converts the spectra array to a MATLAB table.
- plotbygroup, ploterrbygroup - plot spectra by groups using metadata
- plots, ploterrs - plot spectra in new figure with figure formatting options
- peakdecomp - peak decomposition of spectra (formerly gaussdecomp)

Updated methods

- plot can now dynamically generate legends using the LegendText and LegendFun arguments
- proptable (pt) returns custom metadata as well as built-in properties
- sum and mean can now accept a dimension argument as the original MATLAB functions
- find, findindex can search in custom metadata as well

Depreciated methods

- addmeta, removemeta
- gaussdecomp

v. 2.3

New methods
- ploterror plots spectra with errors as shaded areas
- fwhm, to calculate full width at half maximum
- saveh5 for saving data in HDF5 file format
- gaussdecomp for global gaussian decomposition analysis of spectra
- islocalmax and islocalmin to find maxima/minima

Updated methods
- plot returns a handle to the line objects it creates
- splitop - split-apply operations can now use the categorical index within
          - filter argument removed (apply filters before calling splitof)
          - can now call functions with no output arguments
- setx updated to allow for different interpolation methods (for example spline)

v. 2.2

- New methods trapz and cumtrapz for trapezoid integration of spectra
- New method fitbaseline for polynomial fit of spectra / baselines
- The Spectromatic App uses MATLAB 2019b runtime allowing figure export and copying
