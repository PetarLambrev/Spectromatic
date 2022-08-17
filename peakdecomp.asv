function [fitCoeff,linCoeff,fitSpectrum,resid,fitComponents] = peakdecomp(Spectrum,Start,Lower,Upper,peakFunction,nlFitOptions,linFitOptions)
% GAUSSDECOMP Peak Decomposition of spectra
%
% Synthax 
%   [fitCoeff,linCoeff,fitSpectrum,resid,fitComponents] = peakdecomp(Spectrum,Start,Lower,Upper,peakFunction)
%   [fitCoeff,linCoeff,fitSpectrum,resid,fitComponents] = peakdecomp(Spectrum,Start,Lower,Upper,peakFunction,nlFitOptions,linFitOptions)
%
% Description
%   Fit n peaks to a spectrum. The peaks can be modelled as gaussian, skew gaussian
%   or custom function.
%
%
% Input arguments
%   Spectrum - specdata or struct with X and Y fields
%   Start - parameter starting values as an [n,m] matrix where n is the
%      number of components, and m number of parameters (position, width)
%   Lower - parameter lower bounds, must have the same size as Start
%   Upper - parameter upper bounds, must have the same size as Start%  
%   peakFunction - 'gauss','skewgauss', or anonymous function
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

ncomp = height(Start);
npar = width(Start);

if ~exist("nlFitOptions","var") || ~isempty(nlFitOptions)
    nlFitOptions = optimoptions(@lsqnonlin);
    nlFitOptions.MaxFunctionEvaluations = 1000;
    nlFitOptions.MaxIterations = 100;
    nlFitOptions.FunctionTolerance = 1e-6;
end

X = Spectrum.X;

fitCoeff = lsqnonlin(@iterate,Start,Lower,Upper,nlFitOptions);
[resid,linCoeff,fitCurve,fitComponents] = iterate(fitCoeff);
fitCoeff = reshape(fitCoeff,npar,ncomp)';

fitSpectrum = Spectrum;
fitSpectrum.Y = fitCurve;
if isa(Spectrum,"specdata")
    fitSpectrum.ID = string(Spectrum.ID)+" fit";
    fitSpectrum.History = string(Spectrum.History)+"; peakdecomp fit";
end

%% nested functions

    function sim = calculate(allpar)        
        sim = zeros(numel(X),ncomp);        
        for k = 1:ncomp
            p = allpar(k,:);
            
            if isstring(peakFunction) || ischar(peakFunction)
            switch peakFunction
                case "gauss"
                    sim(:,k) = exp(-(X-p(1)).^2/(2*p(2)));
                case "skewgauss"
                    x = (X-p(1))/p(2);
                    sim(:,k) = exp(-x.^2/2);
                    sim(:,k) = sim(:,k) .* (1 + erf(p(3)*x/2));
                otherwise
                    error('%s is not a valid function name',peakFunction)
            end
            else
                sim(:,k) = peakFunction(p,X);
            end
        end
        
    end

    function [res,lincoeff,fitcurve,sim] = iterate(p)        
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
        
        nspec = numel(Spectrum);
        lincoeff = zeros(nspec,ncomp);
        fitcurve = zeros(numel(X),nspec);
        res = zeros(numel(X),nspec);
        
        for k = 1:numel(Spectrum)
            Y = Spectrum(k).Y;
            lincoeff(k,:) = lsqlin(sim,Y,[],[],[],[],LLB,LUB,[],lsqopt);
            fitcurve(:,k) = sim*lincoeff(k,:)';
            res(:,k) = Y - fitcurve(:,k);
        end
    end
end