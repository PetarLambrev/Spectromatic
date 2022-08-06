function B = catfind( A, valueset, catnames)
%CATFIND Find categories in strings
%
% Syntax
%    B = catfind( A, valueset)
%    B = catfind( A, valueset, catnames)
%
% Description
%    B = catfind( A, valueset) searches for keywords in the array A and
%    returns a categorical array of matches. The keywords are specified in 
%    valueset.
%
%    If an element in A matches more than one keyword in valueset, the
%    last (alphabetically) matching keyword will be set as a category of that
%    element.
%    

valueset = string(valueset);
if nargin>2
    catnames = string(catnames);
end
C = repmat("",size(A));
valuesorted = sort(valueset);
for k = 1:numel(valueset)
    cmatch = contains(A, valuesorted(k));    
    if any(cmatch)
        C(cmatch) = valuesorted(k);
    end
end
if nargin>2
    B = categorical(C, valueset, catnames);
else
    B = categorical(C, valueset);
end

