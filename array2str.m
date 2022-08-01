function strarray = array2str(array, delim)
% cellcat concatenate array columns into single strings
%
% strarray = array2str(array)
% strarray = array2str(array, delim)
%
% Concatenates elements in each row in array into a single string.
% using delimiter delim (default is space)
%
% If array is a cell array with only one row
% then each cell must be an array with equal number of elements.

if ~exist('delim','var')
    delim = ' ';
end

n = size(array,1);

if n == 1  && ~ischar(array{1})
    n = numel(array{1});
    newarr = cell(n, size(array,2));
    for l = 1:size(newarr,2)
        if isnumeric(array{l})
            c = num2cellstr(array{l});
        else
            c = cellstr(array{l});
        end
        newarr(:,l) = c(:);
    end
    array = newarr;
end

strarray = cell(n,0);

for k = 1:n
    strarray{k,1} = char(array(k,1));
    for l = 2:size(array,2) 
        strarray{k,1} = [strarray{k}, delim, char(array(k,l))];
    end
end
