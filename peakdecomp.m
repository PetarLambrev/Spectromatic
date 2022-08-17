function [fitCoeff,linCoeff,fitSpectrum,resid,fitComponents] = peakdecomp(Spectrum,Start,Lower,Upper,peakFunction,nlFitOptions,linFitOptions)
% PEAKDECOMP Peak Decomposition of spectra
%
% Synthax 
%   [fitCoeff,linCoeff,fitSpectrum,resid,fitComponents] = peakdecomp(Spectrum,Start,Lower,Upper,peakFunction)
%   [fitCoeff,linCoeff,fitSpectrum,resid,fitComponents] = peakdecomp(Spectrum,Start,Lower,Upper,peakFunction,nlFitOptions,linFitOptions)
%
% Description
%   Fit peaks to a spectrum. The peaks can be modelled as gaussian, skew
%   gaussian, or custom function. Uses nonlinear least squares fit to find
%   the peak shape parameters and linear fit for the peak amplitudes.
%
%
% Example
%   A = specdata.load('IRspectrum_01.txt');
%   S = [1600 100; 1650 100; 1700 100];
%   L = [1550  80; 1600  80; 1650  80];
%   U = [1650 200; 1700 200; 1750 200];
%   [fitCoeff,linCoeff,fitSpec] = peakdecomp(A,S,L,U,'gauss');
%
% Input arguments
%   Spectrum - specdata or struct with X and Y fields
%   Start - parameter starting values as an [n,m] matrix where n is the
%      number of components, and m number of parameters (position, width)
%   Lower - parameter lower bounds, must have the same size as Start
%   Upper - parameter upper bounds, must have the same size as Start%  
%   peakFunction - 'gauss', 'skewgauss', 'lorentz' or anonymous function
%   nlFitOptions - nonlinear fit options object or struct
%   linFitOptions - linear fit options struct
%
% Output
%   fitCoeff - nonlinear fit parameters [n,m] matrix
%   linCoeff - linear coefficients (amplitudes) - n-element size array
%   fitSpectrum - specdata object where Y is the fit curve
%   resid - residuals
%   fitComponents - matrix with column corresponding to the individual
%   peak spectra (unscaled)

% Initial parameters
if numel(Spectrum) > 1
    Spectrum = Spectrum(1);
    warning('Only the first spectrum in the array is analyzed.')
end

ncomp = height(Start); % number of peaks
npar = width(Start);   % number of peak shape parameters

if ~exist("nlFitOptions","var") || ~isempty(nlFitOptions)
    nlFitOptions = optimoptions(@lsqnonlin);
    nlFitOptions.MaxFunctionEvaluations = 1000;
    nlFitOptions.MaxIterations = 100;
    nlFitOptions.FunctionTolerance = 1e-6;
end

X = Spectrum.X;
Y = Spectrum.Y;

% Run fit
fitCoeff = lsqnonlin(@iterate,Start,Lower,Upper,nlFitOptions);
[resid,linCoeff,fitCurve,fitComponents] = iterate(fitCoeff);
fitCoeff = reshape(fitCoeff,npar,ncomp)';

% Create fit spectrum
fitSpectrum = Spectrum;
fitSpectrum.Y = fitCurve;
if isa(Spectrum,"specdata")
    fitSpectrum.ID = string(Spectrum.ID)+" fit";
    fitSpectrum.History = string(Spectrum.History)+"; peakdecomp fit";
end

%% nested functions
    function sim = calculate(allpar)        
        % Objective function
        sim = zeros(numel(X),ncomp);        
        for k = 1:ncomp
            p = allpar(k,:);
            
            if isstring(peakFunction) || ischar(peakFunction)
            switch lower(peakFunction)
                case "gauss"
                    y = exp(-(X-p(1)).^2/(2*p(2)));
                case "skewgauss"
                    x = (X-p(1))/p(2);
                    y = exp(-x.^2/2);
                    y = y .* (1 + erf(p(3)*x/2));
                case "lorentz"
                    y = p(2)^2 ./ ((X-p(1)).^2+p(2)^2);
                otherwise
                    error('%s is not a valid function name',peakFunction)
            end
            else
                y = peakFunction(p,X);
            end
            sim(:,k) = y;
        end
        
    end

    function [resid,linCoeff,fitCurve,sim] = iterate(p)
        % Iteration function
        sim = calculate(p);
        
        if exist("linFitOptions","var") && isfield(linFitOptions,'LLB')
            LLB = linFitOptions.LLB;
        else
            LLB = zeros(1,ncomp);
        end
        if exist("linFitOptions","var") && isfield(linFitOptions,'LUB')
            LUB = linFitOptions.LUB;
        else
            LUB = [];
        end
        lsqopt = struct; lsqopt.Display = 'off';
        linCoeff = lsqlin(sim,Y,[],[],[],[],LLB,LUB,[],lsqopt)';
        fitCurve = sim*linCoeff';
        resid = Y - fitCurve;
    end
end