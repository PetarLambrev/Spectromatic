Spectr-O-Matic Toolbox for MATLAB(R)

by Petar H. Lambrev (2011-2022)

Version 2.3

Spectr-O-Matic is a MATLAB toolbox for analysis of spectroscopy data.  

Requires MATLAB R2017b or later.

For introduction and tutorials, see
https://plambrev.wixsite.com/spectromatic

Download the latest version at
http://www.mathworks.com/matlabcentral/fileexchange/32828-spectr-o-matic

Release notes

v. 2.4

New methods

- table converts the spectra array to a MATLAB table.

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
