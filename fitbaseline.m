function [fitres, gof]  = fitbaseline(SP, n, mask, varargin)
% FITBASELINE - fit polynomial function to spectrum excluding masked regions
%
% Synthax
%    fitres = fitbaseline(SP, n)
%    fitres = fitbaseline(SP, n, mask)
%    fitres = fitbaseline(SP, n, mask, FitOptions)
%    [fitres, gof] = fitbaseline(...)
%
% Description
%   fitres = fitbaseline(SP, n) fits n-th order polynomial to the spectrum
%      in SP.
%
%   fitres = fitbaseline(SP, n, mask) - excludes the masked regions  
%      mask is a m * 2 matrix with each row specifying the starting and 
%      ending X value of a masked region.
%
%
% Example
%   fitres = fitbaseline(SP, 2, [400 450; 600 650])
%      Fits a 2nd order polynomial to the spectrum in SP excluding the 
%      wavelength regions X = 400-450 and X = 600-650. The output fitres
%      is a cfit fit object.
%
% see also fit, cfit

%% Initialize
    FitType = fittype(['poly', num2str(n)]);
    FitOptions = fitoptions(FitType);
     
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
if nargin > 3
    [fitres, gof] = fit(SP.X, SP.Y, FitType, FitOptions, varargin{:});
else
    [fitres, gof] = fit(SP.X, SP.Y, FitType, FitOptions);
end
