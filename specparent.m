classdef specparent
    % specparent - Common methods for 1D and 2D data
    % Spectr-O-Matic version 2.4
    %
    % Parent class for specdata and specdata2D
    %
    % Petar Lambrev, 2012-2022
        
    properties
        ID string = "";        % string identifying the spectrum
        DateTime datetime = datetime;  % date and time
        ExpID string = "";     % experiment ID (default ? current directory)
        X double = [];         % vector (array) of X values
        Y double = [];         % vector (array) of Y values
        XType string = "X";    % X-axis quantity
        XUnit string = "";     % X-axis unit
        YType string = "Y";    % Y-axis quantity
        YUnit string = "";     % Y-axis unit
        Comment string = "";   % various comments
        History string = "";   % operations log
        Metadata struct = struct.empty; % metadata
    end
    
    methods
       %% Queries %%
       function T = metatable(SP,Vars)
           % METATABLE Get all Metadata as a table
           %
           % Synthax
           %   T = metatable(SP)
           %   T = metatable(SP,Vars)
           %   T = mt(...)
           %
           % Description
           %  metatable is similar to proptable but only returns custom
           %  metadata. It may be faster than proptable.
           %
           % See also: PROPTABLE     
           
           if numel(SP) > 1

               % Collect meta variables from ALL spectra
               MetaVars = {};
               for k = 1:numel(SP)
                   MetaVars = union(MetaVars,fieldnames(SP(k).Metadata));
               end

               % Add missing variables and assign as <undefined>
               for k = 1:numel(SP)
                   M = setdiff(MetaVars,fieldnames(SP(k).Metadata));
                   if ~isempty(M)
                       for l = 1:numel(M)
                           SP(k) = SP(k).addmd(M{l},categorical({''}));
                       end
                   end
               end
           end

           MetaStruct = [SP.Metadata];
           if isempty(MetaStruct)
               T = table.empty;
           else
               T = struct2table(MetaStruct,'AsArray',true);
               IDs = [SP.ID];
               IDs(IDs=="") = "noID";
               T.Properties.RowNames = matlab.lang.makeUniqueStrings(IDs);
               if exist('Vars','var')
                   if ~iscell(Vars) && ~isstring(Vars)
                       error('Vars must be a cell array or a string');
                   end
                   T = T(:,Vars);
               end
           end
       end

       function T = mt(SP,VarNames)
           % METATABLE Get all Metadata as a table
           %
           % Synthax
           %   T = metatable(SP)
           %   T = metatable(SP,VarNames)
           %   T = mt(...)
           %
           % Description
           %  metatable is similar to proptable but only returns custom
           %  metadata. It may be faster than proptable esp. for big data.
           %
           % See also: PROPTABLE
           
           if nargin < 2
               T = metatable(SP);
           else
               T = metatable(SP,VarNames);
           end
       end
       
       function T = table(SP,Options)
           % TABLE Convert specdata to table
           %
           % Synthax
           %    T = table(S)
           %    T = table(S,'ExpandMetadata',true)
           %
           % Description
           % 
           % T = table(S)
           % Creates a table with rows for every spectrum and columns
           % containing properties and metadata, such as ID, XType, etc.
           % The actual data (X, Y) are also included in the table.
           %
           % T = table(S,'ExpandMetadata',true)
           % Creates separate columns for custom metadata properties
           arguments
               SP specparent {mustBeNonempty}
               Options.ExpandMetadata logical = false
           end
           
           props = fieldnames(SP);
           
           % Create a blank table
           nprop = numel(props); % number of columns (properties)
           nspec = numel(SP); % number of rows (spectra)
           propcell = cell(nspec,nprop);
           
           for k = 1:numel(SP)
               for p = 1:nprop
                   propcell{k, p} = SP(k).(props{p});
               end
           end
           
           T = cell2table(propcell, 'VariableNames', props);   
           
           % Expand Metadata
           if Options.ExpandMetadata
               T = removevars(T,'Metadata');
               MetaData = metatable(SP);
               T = [T, MetaData];
           end
       end
       
       function T = proptable(SP, PropNames)
           % PROPTABLE Table of properties
           %
           % Synthax
           %    T = proptable(SP)
           %    T = SP.pt
           %    T = SP.pt(props)
           %
           % Description
           % Creates a table with rows for every spectrum and columns
           % containing properties (metadata), such as ID, XType, etc.
           % 
           % T = SP.pt returns ID, DateTime and ExpID only
           % T = SP.pt(props) returns the specified properties
           %
           % Example
           %
           % props = proptable(dat)
           %
           % props = 
           %
           %           ID                  DateTime                  ExpID           dim       XType        XUnit    YType    YUnit        Comment               History         
           %     _______________    ______________________    ___________________    ___    ____________    _____    _____    ______    _____________    ________________________
           % 
           %     'CD ID144.1-1'     ' 2017-03-09 11:22:09'    'D:\...\2017-03-09'    801    'Wavelength'    'nm'     'CD'     'mdeg'    'CD ID141.1'     'load CD ID144.1-1.txt' 
           %     'CD ID144.2-1'     ' 2017-03-09 11:27:44'    'D:\...\2017-03-09'    801    'Wavelength'    'nm'     'CD'     'mdeg'    'CD ID141.2'     'load CD ID144.2-1.txt' 
           %     'CD ID144.3-1'     ' 2017-03-09 11:33:14'    'D:\...\2017-03-09'    801    'Wavelength'    'nm'     'CD'     'mdeg'    'CD ID141.3'     'load CD ID144.3-1.txt' 
           %     'CD baseline-1'    ' 2017-03-09 11:38:37'    'D:\...\2017-03-09'    801    'Wavelength'    'nm'     'CD'     'mdeg'    'CD baseline'    'load CD baseline-1.txt'
           
           IgnoreProps = {'dim','X','Y','Z','Comment','History','Metadata'};
           MetaData = metatable(SP);
           AllProps = table(SP);
           if nargin < 2
               PropNames = fieldnames(SP);                              
               KeepProps = ~matches(PropNames,IgnoreProps);
               PropNames = PropNames(KeepProps);                              
               T = [AllProps(:,PropNames), MetaData];
           else               
               T = [AllProps, MetaData];
               T = T(:,PropNames);           
           end
       end
       
       function T = pt(SP,VarNames)
           % See also: proptable
           if nargin<2
               T = proptable(SP);
           else
               T = proptable(SP,VarNames);
           end
       end
       
       function res = xind(SP,x)
           % XIND index of X values
           %
           % Synthax
           %    xi = A.xind(x)
           %
           % Returns the index of the X value(s) specified by x. 
           % If no exact match is found, returns the index of the closest X
           % If A is an array, xi will be a matrix
           % with rows corresponding to spectra and columns for each x
           res = zeros(numel(SP),numel(x));
           for k = 1:numel(SP)
               for l = 1:numel(x)
                   dif = abs(SP(k).X-x(l));
                   [~, res(k,l)] = min(dif(:));
               end
           end %for
       end %xind
              
       function [maxtab, mintab] = peaks(SP,delta)
           % PEAKS find local extrema
           %
           % Synthax
           %    [maxima, minima] = MySpectrum.peaks(delta)
           % 
           % Description
           % Returns the local minima and maxima 
           % with local amplitude larger than delta.
           % Uses the function peakdet by Eli Billauer.
           %
           % See also: markpeaks

           [maxtab, mintab] = peakdet(SP.Y,delta,SP.X);
       end %peakdet
       
       function res = Yx(SP,Xval)
           % YX Y value at X = Xval
           %
           % Synthax
           %    Y = MySpectrum.Yx(Xval)
           %
           % Description
           % If MySpectrum is an array, Y will also be an array 
           % containing Y values for each spectrum in MySpectrum.
           %
           % This function is useful for calculating spectral parameters,
           % such as peak ratios/differences, etc., for many spectra at once. 
           % 
           % Example
           % Rel = (Data.Yx(X1) ? Data.Yx(X2))./ Data.Yx(X2)
           %
           % calculates the relative ratio [y(x1) - y(x2)] / y(x2)           
           Yn = size(SP(1).Y,2);
           res = zeros(numel(SP),Yn);
           for i = 1:numel(SP)
               if (size(SP(i).Y,2) ~= Yn)
                   error('The spectra in the collection have different dimensions');
               end
               [~,xind] = min(abs(SP(i).X - Xval));
               if (isempty(xind)), error('X value not found in spectrum'); end
               res(i,:) = SP(i).Y(xind,:);
           end %for
       end %Yx       
       
       function res = yatx(SP,Xval)
           % see Yx
           res = Yx(SP,Xval);
       end
       
       %% Arithmetic operators
       function res = plus(SP1,SP2)
           % PLUS Add spectra
           % 
           % Synthax
           % Result = MySpectrum1 + MySpectrum2
           %    MySpectrum1 can be a collection of spectra.
           %    MySpectrum2 can be a single spectrum or
           %    a collection of the same number of spectra as MySpectrum1
           %    or the number of spectra in 1 must be a multiply of the
           %    number of spectra in 2.           
           %
           % Result = MySpectrum1 + Const
           %    Adds constant(s) to a spectrum or array of spectra.
           %    Const must be a single number or a vector with size
           %    equal to the number of spectra in MySpectrum1
           %
           % See also: minus, times, rdivide           
           
           res = binary_arithmetic(SP1,SP2,@plus);
           
       end %plus
       
       function res = minus(SP1,SP2)
           % MINUS Subtract spectra
           % 
           % Synthax
           % Result = MySpectrum1 - MySpectrum2
           %    Subtracts two spectra or arrays of spectra.
           %    MySpectrum2 must contain either a single spectrum or
           %    the same number of spectra as MySpectrum1
           %    or the number of spectra in 1 must be a multiply of the
           %    number of spectra in 2.           
           %
           % Result = MySpectrum1 - Const
           %    Subtracts constant(s) from a spectrum or array of spectra.
           %    Const must be a single number or a vector with size
           %    equal to the number of spectra in MySpectrum1
           %
           % See also: uminus, plus, times, rdivide           
           res = binary_arithmetic(SP1,SP2,@minus);
       end %minus
       
       function res = times(SP1,SP2)
           % TIMES Multiply spectra
           % 
           % Synthax
           % Result = MySpectrum1 * MySpectrum2
           %    Multiply two spectra or arrays of spectra.
           %    MySpectrum2 must contain either a single spectrum or
           %    the same number of spectra as MySpectrum1
           %    or the number of spectra in 1 must be a multiply of the
           %    number of spectra in 2.           
           %
           % Result = MySpectrum1 * Const
           %    Multiplies spectra with constant(s).
           %    Const must be a single number or a vector with size
           %    equal to the number of spectra in MySpectrum1
           %
           % See also: minus, plus, rdivide           
           res = binary_arithmetic(SP1,SP2,@times);
       end %times
       
       function res = log(SP)
           % LOG Logarithm of spectra
           % 
           % Synthax
           % Result = log(MySpectrum)
           %    Convert Y to natural logarithm
    
           res = SP;
           for i = 1:numel(SP)
               res(i).Y = log(SP(i).Y);
               res(i).History = sprintf('log %s',SP(i).History);
           end

       end %log
       
       function res = log10(SP)
           % LOG10 Logarithm of spectra
           % 
           % Synthax
           % Result = log10(MySpectrum)
           %    Calculate the logarithmic Y value.
    
           res = SP;
           for i = 1:numel(SP)
               res(i).Y = log10(SP(i).Y);
               res(i).History = sprintf('log10 %s',SP(i).History);
           end

       end %log
              
       function res = mtimes(SP1,SP2)
           % See also: times
           res = binary_arithmetic(SP1,SP2,@times);
       end %mtimes

       function res = power(SP,P)
           % POWER Power of spectra
           % 
           % Synthax
           % Result = MySpectrum ^ n
           %    Calculate the n-th power of Y.
           if ~isscalar(P)
               error('the second argument must be a scalar value');
           end
    
           res = SP;
           for i = 1:numel(SP)
               res(i).Y = SP(i).Y .^ P;
               res(i).History = sprintf('%s ^ %g',SP(i).History, P);
           end
       end
           
       function res = mpower(SP,P)
           % See also: power
           res = power(SP,P);
       end       
       
       function res = rdivide(SP1,SP2)
           % RDIVIDE Divide spectra
           % 
           % Synthax
           % Result = MySpectrum1 / MySpectrum2
           %    Divides two spectra or arrays of spectra.
           %    MySpectrum2 must contain either a single spectrum or
           %    the same number of spectra as MySpectrum1.
           %    or the number of spectra in 1 must be a multiply of the
           %    number of spectra in 2.
           %
           % Result = MySpectrum1 / Const
           %    Divides constant(s) from a spectrum or array of spectra.
           %    Const must be a single number or a vector with size
           %    equal to the number of spectra in MySpectrum1
           %    or the number of spectra must be a multiply of the number
           %    of constants.
           
           % See also: minus, plus, times
           
           res = binary_arithmetic(SP1,SP2,@rdivide);
       end %rdivide

       function res = mrdivide(SP1,SP2)
           % See also: rdivide
           res = binary_arithmetic(SP1,SP2,@rdivide);
       end %mrdivide
       
       function res = uminus(SP)
           % Unary minus
           % 
           % Synthax
           %    Result = -MySpectrum1 
           %
           % Inverses the sign of the Y values in MySpectrum1
           %
           % See also: minus, plus, times, rdivide           
           res = binary_arithmetic(SP,-1,@times);
       end %uminus
       
       function res = sum(SP,dim)
           % SUM Sum of spectra
           % 
           % Synthax
           %    res = sum(SP)
           %    res = sum(SP,dim)
           %    res = sum(SP,'all')
           %
           % Description
           %   res = sum(SP) sums all spectra across the first
           %   non-singleton dimension
           %   
           %   res = sum(SP,dim) sums spectra across dimension dim
           %
           %   res = sum(SP,'all') sums spectra across all dimensions
           %
           % See also: plus, mean
           if ~exist('dim','var')
               dim = [];
           end           
           res = array_fun(SP,@sum,dim);
       end %sum
       
       function res = mean(SP,dim)
           % MEAN Average spectra
           %
           % Synthax
           %    res = mean(SP)
           %    res = mean(SP,dim)
           %    res = mean(SP,'all')
           %
           % See also: sum
           if ~exist('dim','var')
               dim = [];
           end           
           res = array_fun(SP,@mean,dim);
       end %mean
       
       function [maxY, maxX] = max(SP,Xlim)           
           % MAX Global maxima of spectra
           %
           % Synthax
           %   Y = max(A)
           %   Y = A.max
           %   Y = A.max(XRange)
           %   [Y, X] = A.max(...)
           %
           % Description
           %   Y = max(A) or Y = A.max returns the maximal magnitude 
           %   (Y value) of the spectrum A. 
           %   If A is an array of spectra, Y is a column vector 
           %   of magnitudes for every spectrum in A.
           %   If A is an array of specdata2D, Y is a matrix.
           %
           %   Y = A.max(XRange) finds the global maximum in a range of
           %   X values XRange = [Start, End].
           %
           %   [Y, X] = A.max(...) returns the magnitudes and positions 
           %   (X values) of the maxima.
                      
           Yn = size(SP(1).Y,2);
           maxY = zeros(numel(SP),Yn);
           ix   = zeros(numel(SP),Yn);
           maxX = zeros(numel(SP),Yn);
           for i = 1:numel(SP)
               if (size(SP(i).Y,2) ~= Yn)
                   error('The spectra in the collection have different dimensions');
               end               
               if exist('Xlim','var')
                  SP(i) = SP(i).setxlim(Xlim);
               end
               if numel(SP(i).Y) < 1
                   error('No values found in specified X range. ID = %s',SP(i).ID)
               end
               [maxY(i,:), ix(i,:)] = max(SP(i).Y);
               for j = 1:size(SP(i).Y,2)
                   maxX(i,j) = SP(i).X(ix(i,j));
               end
           end
       end %max
           
       function [minY, minX] = min(SP,Xlim)
           % MIN Global minima of spectra
           %
           % Synthax
           %   Y = min(A)
           %   Y = A.min
           %   Y = A.min(XRange)
           %   [Y, X] = A.min(...)
           %
           % Description
           %   Y = min(A) or Y = A.min returns the minimal magnitude 
           %   (Y value) of the spectrum A. 
           %   If A is an array of spectra, Y is a column vector 
           %   of magnitudes for every spectrum in A.
           %   If A is an array of specdata2D, Y is a matrix.
           %
           %   Y = A.min(XRange) finds the global minimum in a range of
           %   X values XRange = [Start, End].
           %
           %   [Y, X] = A.min(...) returns the magnitudes and positions 
           %   (X values) of the minima.    
           
           Yn = size(SP(1).Y,2);
           minY = zeros(numel(SP),Yn);
           ix   = zeros(numel(SP),Yn);
           minX = zeros(numel(SP),Yn);
           for i = 1:numel(SP)
               if (size(SP(i).Y,2) ~= Yn)
                   error('The spectra in the collection have different dimensions');
               end               
               if exist('Xlim','var')
                  SP(i) = SP(i).setxlim(Xlim);
               end
               if numel(SP(i).Y) < 1
                   error('No values found in specified X range. ID = %s',SP(i).ID)
               end
               [minY(i,:), ix(i,:)] = min(SP(i).Y);
               for j = 1:size(SP(i).Y,2)
                   minX(i,j) = SP(i).X(ix(i,j));
               end
           end           
       end %min       
       
       function S = std(A,w,dim)
           % STD Standard deviation of Y
           %
           % Synthax
           %   S = std(A)
           %   S = std(A,w)
           %   S = std(A,w,dim)
           %
           % Description
           %   S = std(A) calculates the standard deviation S of spectra A.
           %   If A is a collection of spectra,
           %   then S is a spectrum whose Y values contain the standard
           %   deviations of the Y values of the spectra in A.
           %
           %   S = std(A,w) uses weighting factors w to calculate S.
           %   S = std(A,w,dim) calculates along dimension dim (1 or 2)
           
           if ~exist('dim','var')
               dim = [];
           end           
           if exist('w','var')
               S = array_fun(A,@std,dim,w);
           else
               S = array_fun(A,@std,dim);
           end

       end %std

       function S = stderr(A,w,dim)
           % STDERR Standard error of Y
           %
           % Synthax
           %   S = stderr(A)
           %   S = stderr(A,w)
           %   S = stderr(A,w,dim)
           %
           % Description
           %   S = stderr(A) calculates the standard error S of spectra A.
           %   If A is a collection of spectra,
           %   then S is a spectrum whose Y values contain the standard
           %   deviations of the Y values of the spectra in A.
           %
           %   S = stderr(A,w) uses weighting factors w to calculate S.
           %   S = stderr(A,w,dim) calculates along dimension dim (1 or 2)
           %
           % see also: STD

           if ~exist('dim','var')
               dim = [];
           end           
           if exist('w','var')
               [S,d] = array_fun(A,@std,dim,w);
           else
               [S,d] = array_fun(A,@std,dim);
           end
           S = S ./ sqrt(size(A,d));
       end

       function res = xmax(SP)
           % XMAX Maximal X value
           res = zeros(size(SP));
           for k = 1:numel(SP)
               res(k,:) = max(SP(k).X);
           end
       end %xmax
           
       function res = xmin(SP)
           % XMIN Minimal X value
           res = zeros(size(SP));
           for k = 1:numel(SP)
               res(k,:) = min(SP(k).X);
           end
       end %xmin     

       function res = abs(SP)
           % abs - absolute values of Y
           res = SP;
           for i = 1:numel(SP)
               res(i).Y = abs(SP(i).Y);
           end
       end %abs       
       
       
       %% Other calculations
       function res = smooth(SP,varargin)
           % SMOOTH smooth spectra
           %
           % Synthax
           %    new = MySpectrum.smooth(options)
           %
           % Description
           % Uses the function smooth from MATLAB curve-fitting toolbox.
           % See the function?s documentation for the available options.
           %
           % Example
           %    new = MySpectrum.smooth('sgolay',12)
          
           res = SP;
           for i = 1:numel(SP)
               for j = 1:size(SP(i).Y,2)
                   res(i).Y(:,j) = smooth(SP(i).Y(:,j),varargin{:});
               end
               res(i).History = sprintf('%s; smooth', SP(i).History);
           end
       end %smooth

       function res = diff(SP,order)
       % DIFF Differentiate spectra
       %
       % Synthax
       %    DerivSpectrum = MySpectrum.diff(order)
       %
       % Description
       % Calculate the derivative of a spectrum (or array of spectra) 
       % of the specified order (default is 1).

           res = SP;
           if(nargin < 2), order = 1; end
           for i = 1:numel(SP)
               for j = 1:size(res(i).Y,2) 
                   for k = 1:order
                       res(i).Y(:,j) = gradient(res(i).Y(:,j));
                   end
               end
               res(i).History = sprintf('%s; diff', SP(i).History);
           end %for           
       end %diff

       function res = int(SP,Xlim)
           % INT Trapezoid integration of spectra
           %
           % see also trapz
           res = trapz(SP,Xlim);
       end
       
       function res = trapz(SP,Xlim)
           % TRAPZ Trapezoid integration of spectra
           %
           % Synthax
           %    Y = SP.trapz(Xlim);
           %
           % Description
           %    SP is an object containing
           %    one or more spectra and Xlim is an optional
           %    parameter setting the X range for integration
           %    Xlim = [Xmin Xmax]
           %
           % see also cumtrapz
           
           Yn = size(SP(1).Y,2);
           res = zeros(numel(SP),Yn);
           for i = 1:numel(SP)
               if (size(SP(i).Y,2) ~= Yn)
                   error('The spectra in the collection have different dimensions');
               end
              if ~exist('Xlim','var')
                   Xlim = [SP(i).X(1) SP(i).X(end)];
               end
               SPtr = SP(i).setxlim(Xlim);
               for j = 1:Yn                   
                   res(i,j) = trapz(SPtr.X,SPtr.Y(:,j));
               end
           end
       end

       function res = cumtrapz(SP,Xlim)
           % CUMTRAPZ Cumulative trapezoid integration of spectra
           %
           % Synthax
           %    Y = SP.cumtrapz(Xlim);
           %
           % Description
           %    SP is an object containing
           %    one or more spectra and Xlim is an optional
           %    parameter setting the X range for integration
           %    Xlim = [Xmin Xmax]
           %
           % see also trapz
           
           Yn = size(SP(1).Y,2);
           res = specdata.empty;
           for i = 1:numel(SP)
               if (size(SP(i).Y,2) ~= Yn)
                   error('The spectra in the collection have different dimensions');
               end
              if ~exist('Xlim','var')
                   Xlim = [SP(i).X(1) SP(i).X(end)];
               end
               SPtr = SP(i).setxlim(Xlim);
               res(i) = SPtr;
               for j = 1:Yn                   
                   res(i).Y(:,j) = cumtrapz(SPtr.X,SPtr.Y(:,j));
               end
           end
       end
       
       function res = bin(SP,step)
           % BIN Box-car aggregation of data points
           %
           % Synthax
           %    NewSpectrum = MySpectrum.bin(step)
           % 
           % Description
           % Bins (aggregates) datapoints using the specified step. 
           % For example if MySpectrum has datapoints at every 0.1 nm,
           % MySpectrum.bin(10) will return a spectrum with datapoints every 1 nm,
           % where each new datapoint is an average of 10 original datapoints.
           % The bin function is useful to reduce the file size and noise 
           % in the data but it also decreases spectral resolution.
           %
           % See also: smooth

           res = SP;           
           for i = 1:numel(SP)
               Yn = size(SP(i).Y,2);
               bins = fix(numel(SP(i).X)/step);
               res(i).X = zeros(bins,1);
               res(i).Y = zeros(bins,Yn);
               for j = 1:bins
                   ia = (j-1)*step + 1;
                   ib = j*step;
                   res(i).X(j) = mean(SP(i).X(ia:ib));
                   for k = 1:Yn
                       res(i).Y(j,k) = mean(SP(i).Y(ia:ib,k));
                   end
               end %for
               res(i).History = sprintf('%s; bin %.f', SP(i).History,step);
           end %for
       end %bin
       
       function res = norm(SP,Xval)           
       % NORM Normalize spectra
       %
       % Synthax
       %    B = A.norm
       %    B = A.norm('min')
       %    B = A.norm(X)
       %    B = A.norm([X1 X2])
       % 
       % Description
       % B = A.norm normalizes the spectra to the global maximum
       % B = A.norm normalizes the spectra to the global minimum
       % B = A.norm(X) divides by the magnitude at position X
       % B = A.norm([X1 X2]) divides by the maximum in the X range [X1 X2]
       
           res = SP;
           if (nargin < 2), Xval = 'MAX'; end
           for i = 1:numel(SP)
               if (ischar(Xval))
                   switch (upper(Xval))
                       case 'MAX'
                           f = max(SP(i).Y);
                       case 'MIN'
                           f = min(SP(i).Y);
                       otherwise
                           error('Operand %s not recognized',Xval);
                   end %switch
               else
                   if isscalar(Xval) 
                       Xi = find(SP(i).X == Xval,1);
                       if (isempty(Xi)), error('X value not found in spectrum.'); end
                       f = SP(i).Y(Xi,:);
                   else
                       f = SP(i).max(Xval);
                   end                       
               end %if
               res(i).Y = SP(i).Y ./ abs(f);
               res(i).History = sprintf('%s; norm(%s)', SP(i).History, num2str(Xval));
               % res(i).YUnit = 'rel.u.';
           end %for
       end %norm

       %% Modify Properties
      function res = setx(SP,newX,method)
           % SETX Convert X axis
           %
           % Synthax
           %    NewSpectrum = MySpectrum.setx(newX)
           %    NewSpectrum = MySpectrum.setx(newX,method)
           %
           % Description
           % NewSpectrum = MySpectrum.setx(newX)
           % Replaces the original X axis values with newX and
           % redistributes the data. The Y values are interpolated 
           % to correspond to the new X array.
           % The setx function is useful for operations on spectra with 
           % different X scaling. To perform such operation, the X axes must be 
           % made identical first, by using an expression of the form
           %    NewSpectrum1 = MySpectrum1.setx(MySpectrum2.X)
           %
           % setx uses linear interpolation to match the Y values 
           % and does NOT extrapolate. If the new X array contains values outside
           % the range of the original X array, the corresponding Y values 
           % will be ?NaN?. To avoid this, make sure the new X array does not 
           % exceed the range of the original X array. If necessary, 
           % use the setxlim method to limit the data range.
           %
           % NewSpectrum = MySpectrum.setx(newX,method)
           %    method (optional) - set the interpolation method ('spline',
           %    'pchip','cubic','nearest'). default is linear.
           %
           % See also: setxlim, interp1
           classname = mfilename('class');
           if(isa(newX,classname))
               newX = newX.X;
           end
           res = SP;
           if ~exist('method','var')
               method = 'pchip';
           end
           for i = 1:numel(SP)
               res(i).X = newX(:);
               res(i).Y = zeros(numel(newX),size(SP(i).Y,2));
               for j = 1:size(SP(i).Y,2)
                   res(i).Y(:,j) = interp1(SP(i).X,SP(i).Y(:,j),newX(:),method);
               end
           end %for
       end %set X

       function B = shiftx(A,shift)
           % SHIFTX Shift X axis
           %
           % Synthax
           %    B = A.shiftx(offset)
           %
           % Description
           % Shifts the spectrum/spectra along the X axis. 
           % If A is a single spectrum, the result is identical to:
           % B = A; B.X = A.X + offset;
           %
           % Input parameters
           % A - specdata object or array
           % offset - magnitude of the shift in X units. Can be a vector
           % with number of elements equal to the number of spectra in A
           %
           % Output
           % B - spectra with shifted X

           B = A;
           if numel(shift)==1 && numel(A) > 1
               shift = repmat(shift, size(A));
           end
           for k = 1:numel(A)
               B(k).X = A(k).X + shift(k);
           end
               
       end %shiftx
       
       function res = setxlim(SP,Xlim)
           % SETXLIM Set minimum and maximum limits for the X values.
           % 
           % Synthax
           %    NewSpectrum = MySpectrum.setxlim(Xlim)
           %
           % Description
           % The parameter Xlim must be an array of two elements: 
           % Xlim = [Xmin Xmax]. 
           % Any datapoints lying outside the specified range are deleted.
           res = SP;
           for i = 1:numel(SP)
               trim = (res(i).X < Xlim(1)) | (res(i).X > Xlim(2)) | isnan(res(i).X);
               res(i).X(trim) = [];
               res(i).Y(trim,:) = [];
           end %for
       end %setxlim
       
       function P = get(SP,propname,arg)
           % GET get property values
           %
           % Synthax
           % P = SP.get(propname)
           % P = SP.get(propname,arg)
           %
           % Description
           % Returns a selected property as a cell array or number array.
           %    P = SP.get('ID');
           %
           % prop = SP.get(propname,arg) gets properties with an input
           %    P = SP.get('Yx',680);
           %
           % The size of P is equal to the size of SP.
           % If SP is a 2-dimensional array, so will be P.

           if exist('arg','var')
               p0 = SP(1).(propname)(arg);
           else
               p0 = SP(1).(propname);
           end
           if isscalar(p0) && isnumeric(p0)               
               P = repmat(p0,size(SP));
           elseif ischar(p0) || isstring(p0)
               P = repmat("",size(SP));
           else
               P = cell(size(SP));
           end
           for k = 1:numel(SP)
               if exist('arg','var')
                   p1 = SP(k).(propname)(arg);
               else
                   p1 = SP(k).(propname);
               end
               if iscell(P)
                   P{k} = p1;
               else
                   P(k) = p1;                   
               end
           end           
           if iscategorical(p0)
               catnames = categories(SP(1).(propname));
               ncats = numel(catnames);
               P = categorical(P,1:ncats,catnames);
           end
       end
       
       function SP = set(SP,varargin)
           % SET set property values 
           %
           % Synthax
           %     NewSpectrum = MySpectrum.set(Property, Values, Property2, Values2, ...)
           %
           % Description
           % Useful for assigning properties to a collection of spectra.
           % Property is the name of the property to be changed.
           % Values is either a single value or a cell array of values.
           % In that case the number of cells corresponds to the number 
           % of spectra in MySpectrum.
           
           args = struct(varargin{:});
           argf = fieldnames(args);
           numEL = numel(args);
           if (numEL>1)
               if (numEL~=numel(SP))
                   error('Incorrect number of parameters');
               end
           end
           for f = 1:numel(argf)
                property = argf{f};
               for i = 1:numel(SP)
                   if (numEL > 1)
                        SP(i).(property) = args(i).(property);
                   elseif numel(args.(property))==numel(SP)
                       SP(i).(property) = args.(property)(i);
                   else
                       SP(i).(property) = args.(property);
                   end
               end
           end
       end
       
       function res = remove_ext(SP)
           % Remove file extension from the ID string
           %
           % Synthax
           %    res = SP.remove_ext;
           % 
           % If the ID of the input SP is 'file1.ext',
           % the output res ID will be 'file1'
           
           OldIDs = SP.get('ID');
           NewIDs = regexprep(OldIDs, '\.(\w{3})$','');
           res = SP.set('ID',NewIDs);
       end
       
       function res = fi(SP,varargin)
           % See findindex
           res = findindex(SP,varargin{:});
       end
              
       %% Search and Index
       function res = findindex(SP,varargin)
            % FINDINDEX Search for keywords and return a logical index
            %
            % Synthax
            %   Indices = SP.findindex(keyword)            
            %   Indices = SP.findindex(property1, keyword1, property2, keyword2, ...)
            %   Indices = SP.findindex(..., 'w')
            %
            % Description
            % Return a logical array
            % of indices specifying the spectra which match the search
            % conditions.
            %
            % To return the actual objects instead of a logical array,
            % use the method 'find'.
            %
            % Input parameters
            % keyword - string to search for in ID
            %            or cell array of strings (combined with OR)
            %
            % property1, property2, ...  ? specdata properties to be searched
            %   e.g. ?XType?, ?ID?, ?ExpID?, ?Comment?, ?Date?
            % keyword1, keyword2, ... ? strings specifying the text to be searched for. 
            % The search is case-sensitive.
            %
            % keyword can be a string or a cell array of strings. If keyword is a
            % cell array of strings, any match found will be returned (logical OR).
            %
            % If more than one property/keyword pair is supplied, they are combined
            % with logical AND, i.e. a match must be found for each term.
            %
            % If there is no match found, the result is an empty array and a warning is issued.
            % 
            % findindex(..., 'w') - match whole words only (by default
            % partial matches are also valid)
            %
            % Example
            %
            % indx = Data.findindex('ID', {'control', 'ctrl'}, 'YType',
            % 'Fluorescence');
            %
            % The command will search the specdata array Data for spectra
            % whose 'YType' property matches 'Fluorescence' AND whose 'ID'
            % matches either 'control' OR 'ctrl'. The object Fcontrol will
            % contain a logical array with 'true' for all matching spectra.            
            % 
            % See also: find, fi
            
            % Get property table
            T = SP.pt;
            PropNames = T.Properties.VariableNames;
            n = numel(SP);
            idxc = true(n,1);
            
            % parse arguments
            switch numel(varargin)
                case 0
                    error('Cannot perform search. No search terms are given.');
                case 1                                                        
                    terms = varargin;
                    num_args = 1;
                    props = cellstr({'ID'});
                    wholewords = 0;
                otherwise
                    if strcmpi(varargin(numel(varargin)),'w')
                        wholewords = 1;
                        num_args = (numel(varargin)-1)/2;
                    else
                        wholewords = 0;
                        num_args = numel(varargin) / 2;
                    end                    
                    % check for even number of args
                    if (num_args ~= fix(num_args))
                        error('Function requires arguments in pairs.');
                    end
                    props = cell(num_args,1); % properties to search in
                    terms = cell(num_args,1); % text to search for
                    
                    for i = 1:num_args
                        props(i) = varargin((i-1)*2+1); % the property to search in
                        if(~ismember(props{i},PropNames)) % check for a valid property name
                            error('%s is not a valid specdata property',props{i});
                        end
                        terms(i) = varargin((i-1)*2+2);                        
                    end
            end
                % perform search
                for i = 1:num_args
                    idx = false(n,1);
                        str = T.(props{i});

                        term = string(terms{i});                        
                           for k = 1:numel(term)    
                               if iscategorical(str)
                                   idx = idx | str==term(k);                                   
                               elseif wholewords
                                   % whole-word search
                                   idx = idx | contains(str, alphanumericBoundary+term(k)+alphanumericBoundary); % whole words only                                   
                               else
                                   % literal search
                                   idx = idx | contains(str, term(k)); % any match
                               end
                           end                   
                           idxc = idxc & idx;
                end
                if ~any(idxc)
                    warning('MATLAB:specdata:NotFound','The search returned zero results.')
                end
                res = idxc;
        end
        
       function res = find(SP,varargin)
            % FIND Search in properties. 
            %
            % Synthax
            %   newSP = SP.find(keyword)
            %   newSP = SP.find(property1, keyword1, property2, keyword2, ...)
            %
            % Description
            % Returns an array of
            % specdata objects that match the search criteria. 
            %
            % To return a logical array of object indices, use the method
            % 'findindex'
            %
            % Input parameters
            %   property1, property2, ...  ? specdata properties to be searched
            %     e.g. ?XType?, ?ID?, ?ExpID?, ?Comment?, ?Date?
            %   keyword1, keyword2, ... ? strings specifying the text to be
            %   searched for.
            %
            % keyword - string to search for in ID
            %            or cell array of strings (combined with OR)
            %
            % property1, property2, ...  ? specdata properties to be searched
            %   e.g. ?XType?, ?ID?, ?ExpID?, ?Comment?, ?Date?
            % keyword1, keyword2, ... ? strings specifying the text to be searched for. 
            % The search is case-sensitive.
            %
            % keyword can be a string or a cell array of strings. If keyword is a
            % cell array of strings, any match found will be returned (logical OR).
            %
            % If more than one property/keyword pair is supplied, they are combined
            % with logical AND, i.e. a match must be found for each term.
            %
            % If there is no match found, the result is an empty array and a warning is issued.
            % 
            % findindex(..., 'w') - match whole words only (by default
            % partial matches are also valid)
            %
            % Example
            %
            % Fcontrol = Data.find('ID', {'control', 'ctrl'}, 'YType',
            % 'Fluorescence');
            %
            % The command will search the specdata array Data for spectra
            % whose 'YType' property matches 'Fluorescence' AND whose 'ID'
            % matches either 'control' OR 'ctrl'. The object Fcontrol will
            % contain all matching spectra.
            %
            % See also: findindex

            resi = findindex(SP,varargin{:});
            if (isempty(resi))
                res = [];
            else
                res = SP(resi);
            end
       end
        
       function index = autoindex(SP,keywords,keynames)
           % AUTOINDEX Automatic index of keywords
           %
           % Synthax
           %
           %   index = data.autoindex(keywords)
           %   index = data.autoindex(keywords, keynames)
           %   index = data.autoindex
           %
           % Description
           % Creates an index structure with logical arrays
           % matching keywords in the ID of spectra
           %
           % The function will create a struct field for each keyword
           % containing a logical array specifying which spectra IDs
           % match the given keyword. The search is case-sensitive.
           %
           % The names of the fields in the index struct can be optionally
           % specified by the cell array keynames. If not specified, field
           % names will be the same as the keywords.
           %
           % If the keywords array is not given, then a list of keywords
           % will be automatically created from all spectra IDs.
           % Keynames must start with a letter and not contain any special
           % characters (only letters and numbers).
           %
           % Input parameters
           % keywords - (optional) a cell array of strings to search for
           % keynames - (optional) field names to assign to keywords
           % 
           % Example:
           %
           % Suppose the object 'data' contains 4 spectra with IDs:
           % 'GrpA Ctrl', 'GrpB Ctrl', 'GrpA Test', 'GroupB Test'.
           % An automatic keywords index can be created with the line
           % 
           % ix = data.autoindex
           %
           % The result ix will be a structure with four fields:
           % 'Ctrl', 'GrpA', 'GrpB', 'Test'.
           % Each field will contain a 1x4 logical array that specifies
           % which spectra match the given keyword. The index can be used
           % to access any desired spectra, for example:
           % 
           % dif = data(ix.Test) - data(ix.Ctrl)
           %
           % will subtract the two 'Ctrl' spectra from the two 'Test'
           % spectra.
           %
           % Indexes can be combined with logical operators, for example
           % data(ix.GrpA | ix.GrpB) will return all spectra and
           % data(ix.GrpA & ix.Ctrl) will specify only one entry.
           %
           % See also: findindex
           
           index = struct;
           % create a keywords list
           SP1 = SP.remove_ext;
           wholewords = false;
           if nargin == 1
               keywords = [];
               wholewords = true;
               for i = 1:numel(SP1)
                   id = SP1(i).ID;
                   % relpace non-alphanumeric chars with spaces
                   % id = regexprep(id,'[^0-9a-zA-Z_\s]',' '); 
                   ids = regexp(id, "[\s\-;]+", 'split');
                   keywords = cat(2,keywords,ids);                   
               end
               keywords = unique(keywords); % sort unique words
               keynames = matlab.lang.makeValidName(keywords);

           elseif nargin == 2
               keynames = matlab.lang.makeValidName(keywords);
           end

           % create an index
           for k = 1:numel(keywords)               
               try  
                   if wholewords
                       idx = SP1.findindex("ID",keywords(k),'w');
                   else                       
                       idx = SP1.findindex("ID",keywords(k));
                   end
                   if any(idx) && ~all(idx)
                       index.(keynames{k}) = idx;
                   end
               catch ex
                if(~strcmp(ex.identifier,'MATLAB:AddField:InvalidFieldName'))
                    rethrow(ex);
                end
               end
           end
       end

       function B = catfind(SP, valueset, catnames)
           %CATFIND Search for keywords and return a categorical array
           %
           % Syntax
           %    B = A.catfind(valueset)
           %    B = A.catfind(valueset, catnames)
           %    B = catfind(A, ...)
           %
           % Description
           %    B = catfind(valueset) searches for keywords in A.ID
           %    and returns a categorical array of matches. 
           %    The keywords are specified in valueset.
           %
           %    If an element in A matches more than one keyword in valueset, the
           %    last (alphabetically) matching keyword will be set as a category of that
           %    element.
           %
           %    B = A.catfind(valueset, catnames) renames the categories
           %    (keywords) according to the list catnames
           %
           %    B = catfind(A, ...) searches for keywords in the array A
           %    A can be a cell array of strings.
           s = string(SP.get('ID'));
           s = reshape(s, size(SP));
        
           if exist('catnames','var')
               B = catfind(s, valueset, catnames);
           else
               B = catfind(s, valueset);
           end               
       end
       
       function cindex = catindex(SP, vars, newvars)
           % CATINDEX Create a categorical variable index table based on spectra IDs
           %
           % Synthax
           %     
           %     cindex = Data.catindex(vars)
           %     cindex = Data.catindex(vars, newval)
           %     cindex = Data.catindex(vardef)
           %
           % Description
           %
           % Categorical variables are dependent variables that define
           % experimental groups, such as 'Gender'. Variables have a 
           % defined set of values, such as 'Male', 'Female'.
           % The struct vars contains fields for each variable. 
           % Each field must be a cell array of strings defining the
           % possible values for the given variable.
           % Alternatively, variables can be defined as a cell array of
           % strings. Each string must contain a semicolon-separated values
           %  
           % Input parameters
           % vars - structure of experimental variables
           % newvars - structure of variables with new set of values
           % a cell array of strings defining variables and values
           %
           % Output
           % cindex - a table with rows for each file in SP and columns
           % (variables) for each variable in vars.
           %
           % Example
           %
           % Let dat is a specdata object containing six spectra
           % with IDs 'groupA_control.dat', 'groupA_treated.dat', ...
           %
           % The following example defines a categorical index with 
           % two variables - SampleGroup and Treatment,
           % and then uses that index with the splitbinop function 
           % to calculate (treatment - control) spectra for each group
           %
           % vars.SampleGroup = { 'groupA', 'groupB', 'groupC' };
           % vars.Treatment = {'control', 'treated'};
           %
           % datx = dat.catindex(vars);
           % disp(datx)
           %                           SampleGroup    Treatment
           %                           ___________    _________
           % 
           %     groupA_control.dat    groupA         control  
           %     groupA_treated.dat    groupA         treated  
           %     groupB_control.dat    groupB         control  
           %     groupB_treated.dat    groupB         treated  
           %     groupC_control.dat    groupC         control  
           %     groupC_treated.dat    groupC         treated             
           %
           % dif = dat.splitbinop(@minus, datx.Treatment == 'treated', ...
           %                              datx.Treatment == 'control', ...
           %                              datx, 'SampleGroup')
           %
           % disp({dif.History}'
           %     'groupA_treated.dat - groupA_control.dat'
           %     'groupB_treated.dat - groupB_control.dat'
           %     'groupC_treated.dat - groupC_control.dat'  
           %
           % To define a default value, which will be assigned if no
           % other values are found, precede it by a question mark:
           %
           %  vars.Variable1 = {'value1', 'value2', '?default'};


           % Initialize 
           bindex = struct;
           n = numel(SP);
           
           % Loop over variables and spectra
           
           % Loop over variables
           for varname = fieldnames(vars)' 
               varname = varname{1};
               if isstring(vars.(varname))
                   vars.(varname) = cellstr(vars.(varname));
               end

               % Check for a default value               
               defaultvar = '';
               for v = 1:numel(vars.(varname))
                   var = vars.(varname){v};
                   if var(1) == '?'
                       var = var(2:end);
                       vars.(varname){v} = var;
                       defaultvar = var;
                   end
               end
               
               % Initialize index
               bindex.(varname) =repmat({defaultvar},n,1);    
               
               % Loop over Spectra
               for k = 1:n 
%                    bindex.(varname){k} = '';
                   varvals = vars.(varname);
                   varvals = sort(varvals(:));
                   if ~iscell(varvals)
                       vars.(varname) = {varvals};
                   end
                   for l = 1:numel(varvals) % loop over variables
                       var = varvals(l);
                       if iscell(var), var = var{1}; end
                       if strfind(upper(SP(k).ID), upper(var))
                           bindex.(varname){k} = var;
                       end
                   end
               end
               bcats = categorical(bindex.(varname));
               catnames = strtrim(vars.(varname));
               catnames(~ismember(catnames,bcats)) = [];
               bindex.(varname) = reordercats(bcats, catnames);
           end
           
           % Replace values with values from newvars
           if exist('newvars','var')
               for varname = fieldnames(vars)'
                   varname = varname{1};
                   if isfield(newvars, varname)
                       if numel(vars.(varname)) ~= numel(newvars.(varname))
                           error('The variables in vars and newvars must have equal number of values.')
                       end
                       if iscell(newvars.(varname)) || isstring(newvars.(varname))
                           % string variables - simply replace
                           % categories                           
                           oldvarcats = vars.(varname);
                           newvarcats = newvars.(varname);
                           ivars = ismember(oldvarcats, unique(bindex.(varname)));                           
                           cindex.(varname) = renamecats(bindex.(varname),oldvarcats(ivars),newvarcats(ivars));
                       elseif isnumeric(newvars.(varname))
                           % numeric variables
                           cindex.(varname) = zeros(numel(bindex.(varname)),1);
                           for k = 1:numel(vars.(varname))
                               cindex.(varname)(bindex.(varname) == vars.(varname){k}) = newvars.(varname)(k);
                           end
                           cindex.(varname)(isundefined(bindex.(varname))) = NaN;
                       end
                   else
                       cindex.(varname) = bindex.(varname);
                   end
               end
           else
               cindex = bindex;
           end
           IDs = string({SP.ID});
           if numel(unique(IDs)) == n
               cindex = struct2table(cindex,'RowNames',IDs);
           else
               cindex = struct2table(cindex);
           end
       end
 
       function res = splitop(SP, fun, cindex, groupvars, varargin)
           % SPLITOP Split spectra into groups and perform operation
           % 
           %    Synthax
           %      res = SP.splitop(fun, cindex)
           %      res = SP.splitop(fun, cindex, groupvars)
           %      res = SP.splitop(fun, cindex, groupvars, ...)
           %
           % Description
           % splitop performs unary operations (e.g. mean) on spectra
           % The spectra in SP will be split into groups based on the
           % categorical variables in cindex. The function 'fun' will be 
           % called separately for each group of matching spectra and the
           % result will be stored in the array res.
           %
           % Use groupvars to select variables to use for grouping. If
           % groupvars is not given, all variables in cindex will be used.
           %
           %
           % Input parameters
           %    SP - container object for all spectra
           %    fun - function handle for the operation to be performed,
           %          e.g. @mean
           %    cindex - categorical index (table)
           %    groupvars  - variables to use from cindex 
           %                (optional, cell array of strings)
           %    ... - optional arguments to pass on to 'fun'
           %
           % Options (use before any function-specific optional arguments):
           % 
           % 'PassGroupVars' :
           % Pass the values of the grouping variables to the function fun
           %
           % 'PassIndices': 
           % Pass a additional logical array specifying which spectra in SP 
           % are passed to the function @fun
           % 
           % 'PassCatIndex':
           % Pass the categorical index (only relevant rows for each group)
           %
           % See also: splitapply, splitbinop, catindex
          
           n = numel(SP);
           if size(cindex, 1) ~= n
               error('The number of rows of cindex must be equal to the number of spectra.')
           end
           if isstruct(cindex)
               cindex = struct2table(cindex);
           end
           vars = string(cindex.Properties.VariableNames);
           if exist('groupvars', 'var') && ~isempty(groupvars)
               groupvars = string(groupvars);
               vars = intersect(vars, groupvars);
           end
           nvar = numel(vars);                                 
           ia = true(n, 1);
           if nargout
               res = {};
           end
               
           while any(ia)
               iaf = find(ia); 
               i1 = iaf(1);
               ix = ia;
               for k = 1:nvar                   
                   if iscategorical(cindex.(vars{k})) 
                       if isundefined(cindex.(vars{k})(i1))                                        
                            ix = ix & isundefined(cindex.(vars{k}));                  
                       else
                            ix = ix & cindex.(vars{k}) == cindex.(vars{k})(i1);
                       end
                   elseif isstring(cindex.(vars{k}))
                       if ismissing(cindex.(vars{k})(i1))                                        
                            ix = ix & ismissing(cindex.(vars{k}));                  
                       else
                            ix = ix & cindex.(vars{k}) == cindex.(vars{k})(i1);
                       end                           
                   elseif ischar(cindex.(vars{k}))
                       ix = ix & strcmp(cindex.(vars{k}), cindex.(vars{k})(i1));
                   elseif isnumeric(cindex.(vars{k})(i1))
                       if isnan(cindex.(vars{k})(i1))
                           ix = ix & isnan(cindex.(vars{k}));
                       else
                           ix = ix & cindex.(vars{k}) == cindex.(vars{k})(i1);
                       end
                   else
                       error('The data type in the table cannot be used for indexing.');
                   end
               end
               ia = ia & ~ix;
               if isempty(varargin)
                   if nargout
                       res = [res; fun(SP(ix))];
                   else
                       fun(SP(ix));
                   end
               else
                   if strcmpi(varargin{1},'PassGroupVars')
                       % Pass grouping variables
                       ngv = numel(groupvars);
                       groupvarlist = cell(1,ngv);
                       for igv = 1:ngv
                           groupvarlist(igv) = cellstr(cindex{find(ix,1,'first'),groupvars(igv)});
                       end
                       if nargout
                           res = [res; fun(SP(ix), groupvarlist)];
                       else 
                           fun(SP(ix), groupvarlist);
                       end
                   elseif strcmpi(varargin{1},'PassIndices')
                       % Pass Index
                       if length(varargin) > 1
                           if nargout
                               res = [res; fun(SP(ix), ix, varargin{2:end})];
                           else
                               fun(SP(ix), ix, varargin{2:end});
                           end
                       else
                           if nargout
                               res = [res; fun(SP(ix), ix)];
                           else
                               fun(SP(ix), ix);
                           end
                       end
                   elseif strcmpi(varargin{1},'PassCatIndex')
                       % Pass categorical index
                       if length(varargin) > 1
                           if nargout
                               res = [res; fun(SP(ix), cindex(ix,:), varargin{2:end})];
                           else
                               fun(SP(ix), cindex(ix,:), varargin{2:end});
                           end
                       else
                           if nargout
                               res = [res; fun(SP(ix), cindex(ix,:))];
                           else
                               fun(SP(ix), cindex(ix,:));
                           end
                       end                       
                   else
                       if nargout
                           res = [res; fun(SP(ix), varargin{:})];
                       else
                           fun(SP(ix), varargin{:});
                       end
                   end
               end
           end                
       end
       
       function res = splitbinop(SP, fun, cindex, group1, group2, groupvars)
           % SPLITBINOP Split spectra into groups and perform binary operation 
           % 
           % Synthax
           %      res = SP.splitbinop(fun, catindex, group1, group2, groupvars)
           %     
           % Input parameters
           %    SP - container object for all spectra
           %    fun - function handle for the operation to be performed,
           %          e.g. @minus
           %    group1 - logical array selecting the first 
           %             group of spectra ([] for all spectra)
           %    group2 - logical array selecting the second
           %             group of spectra 
           %    catindex - categorical index (struct)
           %    groupvars  - variables to use from catindex 
           %                (optional, cell array of strings)
           %
           % Description
           % splitbinop performs binary operations (addition, subtraction, ...)
           % on two arrays of spectra, group1 and group2
           % The function uses the categorical index catindex to find
           % pairs of spectra matching the same category.
           %
           % For each spectrum in group1, the function will search for a
           % matching spectrum in group2 based on the catindex variables
           % specified by groupvars. If groupvars is not given, all
           % variables in catindex are used.
           %
           % If no match is found in group2, the respective spectrum 
           % in group1 is copied without change. 
           % If more than one match is found in group2, the matching
           % spectra are averaged.
           % 
           %
           % See also: splitop, catindex
          
           if isstruct(cindex)
               cindex = struct2table(cindex);
           end
           vars = cindex.Properties.VariableNames;
           if exist('groupvars', 'var')
               vars = intersect(vars, groupvars);
           end
           
           nvar = numel(vars);           
           if isempty(group1), group1 = true(numel(SP),1); end
%            res = [];
           resi = find(group1);
           n = numel(resi);
           
           % convert to categorical
           for v = 1:nvar
               if ~isnumeric(cindex.(vars{v}))
                   cindex.(vars{v}) = categorical(cindex.(vars{v}));
               end
           end
           for k = 1:n
               i1 = group1(:); i2 = group2(:);
               for v = 1:nvar
                   var_value = cindex.(vars{v})(resi(k));
                   if iscategorical(var_value) && isundefined(var_value)
                       i1 = i1 & isundefined(cindex.(vars{v}));
                       i2 = i2 & isundefined(cindex.(vars{v}));
                   elseif ischar(var_value)
                       i1 = i1 & strcmp(cindex.(vars{v}), var_value);
                       i2 = i2 & strcmp(cindex.(vars{v}), var_value);
                   else
                       i1 = i1 & cindex.(vars{v}) == var_value;
                       i2 = i2 & cindex.(vars{v}) == var_value;
                   end
               end
               if any(i1) 
                   if any(i2)
                       res(k) = fun(SP(resi(k)), mean(SP(i2)));                       
                   else
                       warning('Matching spectra not found in group2.');
                   end
               end
           end
       end
     
       function sorted = sort(SP,field)
           % SORT Sort spectra 
           %
           % Synthax
           %   sorted = data.sort(field)
           %   sorted = data.sort
           %
           % Description
           % Sorts the spectra by the specified field (default is ID). 
           %
           %  Example:
           %  data contains IDs of 'ID1', 'ID3', 'ID2'
           %  data.sort will rearrange the spectra as 'ID1', 'ID2', 'ID3'
           
           if nargin < 2
               field = 'ID';
           end
           PropTable = SP.proptable(field);
           [~, index] = sortrows(PropTable);
           sorted = SP(index);
       end

       %% Metadata
       function res = setmetadata(SP,Property,Values)
          % SETMETADATA Assign metatadata to spectra
          %
          % Synthax
          %   res = setmetadata(SP,Property,Values)
          %   res = setmetadata(SP,T)
          %   res = setmd(...)
          %
          % Description
          %   res = setmetadata(SP,T) assigns metadata in table (categorical
          %     index) T to the spectra SP. The number of rows in T must
          %     match the number of spectra in SP.
          %
          %   res = setmetadata(SP,Property,Values) sets a specific
          %     metadata property with the given value(s)

          res = SP;          
          if istable(Property)
              % Set the whole metadata structure
              T = Property;
              if height(T) ~= numel(SP)
                  error("The number of rows in T must match the number of spectra in SP")
              end

              % Implementation
              MetaStruct = table2struct(T);

              for k = 1:numel(SP)
                  res(k).Metadata = MetaStruct(k);
              end

          else
              % Set individual property/values
              if numel(Values) == 1
                  Values = repmat(Values,size(SP));
              end

              for k = 1:numel(SP)                  
                  if isfield(SP(k).Metadata,Property)
                      res(k).Metadata.(Property) = Values(k);
                  else
                      res(k) = res(k).addmd(Property,Values(k));
                  end
              end
          end
       end
       
       function res = setmd(SP,Property,Values)
          % SETMETADATA Assign metatadata to spectra
          %
          % Synthax
          %   res = setmd(SP,T)
          %   res = setmd(SP,Property,Values)
          
          %
          % Description
          %   res = setmd(SP,T) assigns metadata in table (categorical
          %     index) T to the spectra SP. The number of rows in T must
          %     match the number of spectra in SP.
          %
          %   res = setmd(SP,Property,Values) sets a specific metadata
          %     property with the given value(s)
          %
          % SEE ALSO: setmetadata
            
          if exist("Values","var")        
              res = setmetadata(SP,Property,Values);
          else
              res = setmetadata(SP,Property);
          end
    end
          
       function res = metaindex(SP, Vars, newVars)
           % METAINDEX Create categorical index and assign to metadata
           %
           % Synthax
           %   res = metaindex(SP, Vars)
           %   res = metaindex(SP, Vars, newVars)
           %
           % Description
           %   res = metaindex(SP, Vars) creates a categorical index from
           %   the variables in struct Vars and assigns it as Metadata in
           %   the spectra in SP. It is equivalent to the commands:
           %
           %   Idx = catindex(SP,Vars)
           %   res = setmetadata(SP,Idx)
           %
           % See also: CATINDEX, SETMETADATA           
           
           if nargin < 3
               Idx = catindex(SP,Vars);
           else
               Idx = catindex(SP,Vars,newVars);
           end
           res = setmetadata(SP,Idx);           
       end
        
       function res = addmetadata(SP,varargin)
           % ADDMETADATA Add metadata to spectra
           %
           % Synthax           
           %   res = addmetadata(SP,VarName,Values)
           %   res = addmetadata(SP,Var1,Vals1,Var2,Vals2,...)
           %   res = addmetadata(SP,VarStruct)
           %   res = addmd(...)
           %
           % Description
           %   SP = addmetadata(SP,VarName,Values) assigns
           %   Values to the metadata variable VarName in all elements of SP
           %
           %   SP = addmetadata(SP,Var1,Vals1,Var2,Vals2,...) assigns
           %   multiple metadata variables
           %
           %   SP = addmetadata(SP,VarStruct) assigns metadata using the
           %   struct VarStruct containing each variable as a field
           %
           % See also: deletemetadata, metatable, metaindex, find, findindex

           
           if nargin==2 && isstruct(varargin{1})
               varStruct = varargin{1};
           else
               varStruct = struct(varargin{:});
           end
           varNames = fieldnames(varStruct);
           res = SP;
           for varName = varNames'
               for k = 1:numel(SP)
                   varValue = varStruct.(varName{1});
                   res(k).Metadata(1).(varName{1}) = varValue(k);
               end
           end
       end
       
       function res = addmd(SP,varargin)
           % ADDMETADATA Add metadata to spectra
           %
           % Synthax           
           %   res = addmetadata(SP,VarName,Values)
           %   res = addmetadata(SP,Var1,Vals1,Var2,Vals2,...)
           %   res = addmetadata(SP,VarStruct)
           %   res = addmd(...)
           %
           % Description
           %   SP = addmetadata(SP,VarName,Values) assigns
           %   Values to the metadata variable VarName in all elements of SP
           %
           %   SP = addmetadata(SP,Var1,Vals1,Var2,Vals2,...) assigns
           %   multiple metadata variables
           %
           %   SP = addmetadata(SP,VarStruct) assigns metadata using the
           %   struct VarStruct containing each variable as a field
           %
           % See also: deletemetadata, metatable, metaindex, find, findindex      
           
           res = addmetadata(SP,varargin{:});
       end
       
       function res = deletemetadata(SP,VarNames)
           % DELETEMETADATA Delete all metadata from spectra
           %
           % Synthax
           %   res = deletemetadata(SP)
           %   res = deletemetadata(SP,VarNames)
           %   res = deletemd(...)
           %
           % See also: deletemetadata, metatable, metaindex, find, findindex
           
           res = SP;
           for k = 1:numel(res)
               if nargin < 2
                   res(k).Metadata = struct.empty;
               else
                   res(k).Metadata = rmfield(SP(k).Metadata,VarNames);
               end
           end           
       end    
       
       function res = deletemd(SP,VarNames)
           % DELETEMETADATA Delete all metadata from spectra
           %
           % Synthax
           %   res = deletemetadata(SP)
           %   res = deletemetadata(SP,VarNames)
           %   res = deletemd(...)
           %
           % See also: deletemetadata, metatable, metaindex, find, findindex
           
           if nargin < 2
               res = deletemetadata(SP);
           else
               res = deletemetadata(SP,VarNames);
           end
       end
    end %methods
    
    %% Private methods
    methods (Access = protected, Hidden = true)
        function mode = get_arithmetic_mode(SP1, SP2)
           % validate arguments of arithmetic functions
           % SP1 must be collection of spectra
           % SP2 must be either spectra or scalars
           % the number of items in SP1 and SP2 must be equal
           % or SP2 must be singleton
           % or SP1 must have N times the elements of SP2
           % if SP2 are spectra, they must be size-comparable to SP1
           % if the requirements are not fulfilled, generate error
           % otherwise return the type of the second operand:
           %   0: scalar, single
           %   1: scalar, multiple, equal number
           %   2: scalar, multiple, nonequal number
           %   3: object, single
           %   4: object, multiple, equal number
           %   5: object, multiple, nonequal number
           classname = mfilename('class');
           if(isreal(SP1) || isreal(SP2))
               objectmode = 0;
           elseif(isa(SP1,classname) && isa(SP2,classname))
               objectmode = 3;
           else
               error('Unrecognized operand data type. Arithmetic can be done only with objects of the same type or scalars.');
           end %if    
           n1 = numel(SP1); n2 = numel(SP2);
           if(n1==1 || n2==1)
               multimode = 0;
           elseif (n1==n2)                   
               multimode = 1;
           elseif rem(n1,n2) == 0
               multimode = 2;
           else
               error ('Cannot perform arithmetic on objects with different number of spectra.');
           end %if
           mode = objectmode + multimode;
        end
        
        function res = binary_arithmetic(SP1, SP2, fun)
            % perform binary arithmetic
            % SP1, SP2 - operands
            % fun - function handle
            mode = get_arithmetic_mode(SP1, SP2);
            switch func2str(fun)
                case 'plus'
                    op = '+';
                case 'minus'
                    op = '-';
                case 'times'
                    op = '*';
                case 'rdivide'
                    op = '/';
                otherwise
                    op = func2str(fun);
            end
            switch (mode)
                case 0 % SP1/SP2 is scalar, SP2 single
                    if(isreal(SP1))
                        res = SP2;
                        for i = 1:numel(SP2)
                            res(i).Y = fun(SP2(i).Y, SP1);
                            res(i).History = sprintf('%s; %s %g',SP2(i).History,op,SP1);
                        end
                    else   
                        res = SP1;
                        for i = 1:numel(SP1)
                            res(i).Y = fun(SP1(i).Y, SP2);
                            res(i).History = sprintf('%s; %s %g',SP1(i).History,op,SP2);
                        end
                    end
                case 1 % SP1/SP2 is scalar, SP2 multiple
                    if(isreal(SP1))
                        res = SP2;
                        for i = 1:numel(SP2)
                            res(i).Y = fun(SP2(i).Y, SP1(i));
                            res(i).History = sprintf('%s; %s %g',SP2(i).History,op,SP1(i));
                        end
                    else
                        res = SP1;
                        for i = 1:numel(SP1)
                            res(i).Y = fun(SP1(i).Y, SP2(i));
                            res(i).History = sprintf('%s; %s %g',SP1(i).History,op,SP2(i));
                        end
                    end
                    
                case 2 % SP1/SP2 is scalar, SP2 multiple, unequal
                    if numel(SP1) > numel(SP2)
                        res = SP1; opn = SP2;
                    else
                        res = SP2; opn = SP1;
                    end
                    n = numel(opn);
                    for i = 1:n
                        for j = i:n:numel(res)
                            res(j).Y = fun(res(j).Y, opn(i));
                            res(j).History = sprintf('%s; %s %g',res(j).History,op,opn(i));
                        end
                    end
                        
                case 3 % SP1 & SP2 are objects, SP2 single
                    res = SP1;
                    for i = 1:numel(SP1)
                        cmp = comp(SP1(i),SP2);
                        if cmp == 1
                            SP1(i) = SP1(i).setx(SP2);
                            res(i).X = SP1(i).X;
                        elseif cmp == 2
                            error('Cannot perform arithmetic operation with spectra. The spectral axes do not match.');
                        end                        
                        res(i).Y = fun(SP1(i).Y, SP2.Y);
                        res(i).History = sprintf('%s; %s %s',SP1(i).History,op,SP2.ID);
                    end
                    
                case 4 % SP1 & SP2 are objects, SP2 multiple, equal number
                    res = SP1;
                        for i = 1:numel(SP1)
                            cmp = comp(SP1(i),SP2);
                            if cmp == 1
                                SP2(i) = SP2(i).setx(SP1(i));
                                res(i).X = SP2(i).X;
                            elseif cmp == 2
                                error('Cannot perform arithmetic operation with spectra. The spectral axes do not match.');
                            end
                            res(i).Y = fun(SP1(i).Y, SP2(i).Y);
                            res(i).History = sprintf('%s; %s %s',SP1(i).History,op,SP2(i).ID);
                        end
                        
                case 5 % SP1 & SP2 are objects, SP2 multiple, unequal number
                    if size(SP1,2) ~= size(SP2,2)
                        error('The number of spectra or columns of spectra must be equal');
                    else
                        if size(SP2,1) == 1 
                            res = SP1; opn = SP2;
                        elseif size(SP1,1) == 1
                            res = SP2; opn = SP1;
                        else
                            error('The second operand must contain only one row of spectra.')
                        end
                    end
                    [m,n] = size(res);
                    for j = 1:m
                        for i = 1:n
                            cmp = comp(res(j,i),opn(1,i));
                            if cmp == 1
                                res(j,i) = res(j,i).setx(opn(1,i));                                
                            elseif cmp == 2
                                error('Cannot perform arithmetic operation with spectra. The spectral axes do not match.');
                            end
                            res(j,i).Y = fun(res(j,i).Y, opn(1,i).Y);
                            res(j,i).History = sprintf('%s; %s %s',res(j,i).History,op,opn(1,i).ID);
                        end  
                    end
            end % switch
        end        
        
        function [res,dim] = array_fun(SP,fun,dim,weights)
            % AGGREGATEFUN Perform array function over specified dimension            
            SP_size = size(SP);
            if length(SP_size) > 2
                error('Array functions only operate on 1- or 2-dimensional arrays of spectra.');
            end
            if ~exist('weights','var')
                weights = [];
            end
            if isempty(dim)
                if SP_size(1) > 1
                    dim = 1;                    
                elseif SP_size(2) > 1
                    dim = 2;
                else
                    dim = 1;
                end
            elseif ischar(dim) && strcmpi(dim,'all')
                SP = SP(:);
                SP_size = size(SP);
                dim = 1;
            elseif dim < 1 || dim > 2
                error("Dimension can only be 1, 2, or 'all'.")
            end
            
            funInfo = functions(fun);
            differentX = false;
            switch dim
                case 1
                    res = SP(1,:);        
                    for k = 1:SP_size(2)
                        if SP_size(1) > 1
                            x1 = SP(1,k).X;                            
                            for m = 2:SP_size(1)
                                x = intersect(x1,SP(m,k).X);
                                if ~isequal(x,x1)
                                    differentX = true;
                                    x1 = x;
                                end
                            end           
                            if differentX                                
                                SP(:,k) = SP(:,k).setx(x);
                            end
                        end                        
                        Ys = cat(2,SP(:,k).Y);
                        if ismember(funInfo.function,'std')
                            % function with weights
                            res(k).Y = fun(Ys,weights,2);                        
                        else 
                            res(k).Y = fun(Ys,2);                        
                        end
                    end
                case 2
                    res = SP(:,1);
                    for k = 1:SP_size(1)
                        if SP_size(2) > 1
                            x1 = SP(k,1).X;                            
                            for m = 2:SP_size(2)
                                x = intersect(x1,SP(k,m).X);
                                if ~isequal(x,x1)
                                    differentX = true;
                                    x1 = x;
                                end
                            end             
                            if differentX
                                SP(k,:) = SP(k,:).setx(x);
                            end
                        end
                        Ys = cat(2,SP(k,:).Y);
                        
                        if ismember(funInfo.function,'std')
                            % function with weights
                            res(k).Y = fun(Ys,weights,2);                        
                        else 
                            res(k).Y = fun(Ys,2);                        
                        end
                    end 
            end
            if differentX
                warning('X values do not match across spectra. Some data have been discarded.');
            end
        end
            

     end % private methods
    
    %% Abstract methods
    methods (Access = protected, Abstract = true)
        match = comp(SP1, SP2)        
    end

    methods (Static)
        function index = ind(arr, val)
            % ind - Find the index of an array element of a given value
            %
            % index = ind(arr, val)
            %
            % ind searches for the value val in the array arr
            % and returns the index of the closest matching array element
            
            if isempty(arr)
                index = [];
            else
                if isempty(val)
                    error('val cannot be empty.')
                end
                if ~isvector(val)
                    error('val must be a scalar or a vector')
                end
                index = zeros(1, numel(val));
                for k = 1:numel(val)
                    dif = abs(arr-val(k));
                    [~, index(k)] = min(dif(:));
                end
            end
        end    
    end
end
