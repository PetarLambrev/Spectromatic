function fnames = getdir(filter,varargin)
    % Read directory contents
    %
    % filenames = getdir(filepattern)
    % filenames = getdir(filter, Name, Value)
    %
    % returns a cell array of filenames matching the filepattern
    %
    % filepattern can contain path and wildcards
    %
    % Name, Value pairs
    %
    % FullPath - true | false (default) - list file names with full path
    %
    % Example:
    %
    % filenames = getdir('Data\*.txt', 'FullPath', true);
    
    options = struct;
    if nargin > 1
        options = struct(varargin{:});
    end
    if isfield(options, 'FullPath') && options.FullPath
        FullPath = true;
    else 
        FullPath = false;
    end
        
    if FullPath
        [status, files] = fileattrib(filter);
    else
        files = dir(filter);
        status = 1;
    end
    
    if status 
        n = length(files);
    if n>0
       fnames = cell(n,1);
       for i = 1:n
           if FullPath
               fnames(i) = cellstr(files(i).Name);
           else
               fnames(i) = cellstr(files(i).name);
           end
       end
    else
        fnames = {};
    end
    else
        fnames = {};
    end
end