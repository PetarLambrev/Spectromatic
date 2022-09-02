function [B, fitres, gof]  = fitpoly(SP, n, mask, varargin)
% fitpoly - fit polynomial function to spectrum excluding masked regions
%
% Synthax
%    B = fitpoly(SP, n)
%    B = fitpoly(SP, n, mask)
%    B = fitpoly(SP, n, mask, FitOptions)
%    [B, fitres, gof] = fitpoly(...)
%
% Description
%   B = fitpoly(SP, n) fits n-th order polynomial to the spectrum
%      in SP and returns the fit as a specdata object.
%
%   B = fitpoly(SP, n, mask) - excludes the masked regions  
%      mask is a m * 2 matrix with each row specifying the starting and 
%      ending X value of a masked region.
%
%   [B, fitres, gof] = fitpoly(...) returns also a cfit object fitres
%   and gof
%
% Note: before fit the X values will be shifted to zero at their median.
%  The fit baseline B shifts back the baseline to the original X values but
%  the fit coefficients in fitres refer to the centred X values. 
%
% Example
%   fitres = fitpoly(SP, 2, [400 450; 600 650])
%      Fits a 2nd order polynomial to the spectrum in SP excluding the 
%      wavelength regions X = 400-450 and X = 600-650. The output fitres
%      is a cfit fit object.
%
%
% see also fit, cfit

%% Initialize
    FitType = fittype(['poly', num2str(n)]);
    FitOptions = fitoptions(FitType,'Robust','On');
     
    if ~exist('mask','var')
        mask = [];
    end

    if length(SP) > 1
        SP = SP(1);
        warning('Warning. Only first spectrum will be fitted.');
    end
    
    
%% Exclude masked regions
if ~isempty(mask)
    if size(mask,2) ~= 2
        error('The mask must have two columns [start, end].')
    end
    
    FitOptions.Weights = ones(length(SP.X),1);
    for k = 1:size(mask,1)
        Xi = ind(SP.X, mask(k,:));
        FitOptions.Weights(min(Xi):max(Xi))= 0;
    end
end

%% Run fit
Xc = median(SP.X);
if nargin > 3
    [fitres, gof] = fit(SP.X-Xc, SP.Y, FitType, FitOptions, varargin{:});
else
    [fitres, gof] = fit(SP.X-Xc, SP.Y, FitType, FitOptions);
end

B = SP;
B.Y = fitres(B.X-Xc);
B.ID = B.ID + " fit";