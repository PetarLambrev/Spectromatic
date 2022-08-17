function [fitcoeff,lincoeff,fitcurve,resid,comp] = gaussdecomp(Spectra,Start,Lower,Upper,FitOpt)
% GAUSSDECOMP Skew-Gaussian Decomposition of spectra
% 
% Deprecated. Use peakdecomp instead.

% Initial parameters

if isvector(Start)
    ncomp = numel(Start);
else
    ncomp = numel(Start)/3;
end

    NLFitOpt = optimoptions(@lsqnonlin);
    NLFitOpt.MaxFunctionEvaluations = 3000;
    NLFitOpt.MaxIterations = 1000;
    NLFitOpt.FunctionTolerance = 1e-6;

% Convert to wavenumbers
if isvector(Start)    
    xf = 1e7./Start;
    Start = zeros(1,ncomp*3);
    Start(1:3:end) = xf;
    Start(2:3:end) = mean(abs(diff(xf)));
else
    Start(1:3:end) = 1e7./Start(1:3:end);
    if exist('Lower','var')
        Lower(1:3:end)  = 1e7./Lower(1:3:end) ;        
    end
    if exist('Upper','var')
        Upper(1:3:end) = 1e7./Upper(1:3:end);
    end
end
X = 1e7./Spectra(1).X;

fitcoeff = lsqnonlin(@iterate,Start,Lower,Upper,NLFitOpt);
[resid,lincoeff,fitcurve,comp] = iterate(fitcoeff);
fitcoeff = reshape(fitcoeff,3,ncomp)';
fitcoeff(:,1) = 1e7./fitcoeff(:,1);

% sort by wavelength
[~,I] = sort(fitcoeff(:,1));
fitcoeff = fitcoeff(I,:);
lincoeff = lincoeff(I,:);

%% nested functions

    function sim = calculate(p)        
        sim = zeros(numel(X),ncomp);
        for k = 1:ncomp
            m = (k-1)*3;
            x = (X-p(m+1))/p(m+2);
            sim(:,k) = exp(-x.^2/2);
            sim(:,k) = sim(:,k) .* (1 + erf(p(m+3)*x/2));
        end
        
    end

    function [res,lincoeff,fitcurve,sim] = iterate(p)        
        sim = calculate(p);
        
        if isfield(FitOpt,'LLB')
            LLB = FitOpt.LLB;
        else
            LLB = zeros(1,ncomp);
        end
        if isfield(FitOpt,'LUB')
            LUB = FitOpt.LUB;
        else
            LUB = [];
        end
        lsqopt = struct; lsqopt.Display = 'off';
        
        nspec = numel(Spectra);
        lincoeff = zeros(ncomp,nspec);
        fitcurve = zeros(numel(X),nspec);
        res = zeros(numel(X),nspec);
        
        for k = 1:numel(Spectra)
            Y = Spectra(k).Y;
            lincoeff(:,k) = lsqlin(sim,Y,[],[],[],[],LLB,LUB,[],lsqopt);
            fitcurve(:,k) = sim*lincoeff(:,k);
            res(:,k) = Y - fitcurve(:,k);
        end
    end
end