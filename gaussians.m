function Components = gaussians(par, X)
% compute gaussian components
%
% Usage:
%    Components = gaussians(par)
%
% par is either: 
% 1) vector of parameters [a1, b1, c1, a2, b2, c2, ...], or
% 2) cfit object contaning the result of a gaussian model fit
%
% X is a vector of X values for which the function is computed
%
% The function returns a set of SpecData objects calculated as
% Components(1) = a1*exp(-((X-b1)/c1)^2), ...

if isa(par,'cfit')
    coeffs = coeffvalues(par);
else
    if isvector(par)
        if rem(length(par),3)
            error('The argument par must be a vector of the type [a1, b1, c1, ...]');
        end         
    coeffs = par;
    else
        error('The argument par must be a vector of the type [a1, b1, c1, ...]');
    end
end    
ncomp = length(coeffs)/3;
Components = specdata.empty(0,ncomp);

for i = 0:ncomp-1
    a = coeffs(i*3+1);
    b = coeffs(i*3+2);
    c = coeffs(i*3+3);
    Y = zeros(length(X),1);
    for j = 1:length(Y)
        Y(j) = a*exp(-((X(j)-b)/c)^2);
    end
    Components(i+1) = specdata(X,Y,sprintf('g%d',i+1));
end