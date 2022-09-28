classdef specdata < specparent
    % specdata - one-way (X/Y) spectroscopy data
    % Spectr-O-Matic version 2.4
    %
    % Container for X-Y data (wavelength, amplitude)
    % and methods for general spectral manipulation - 
    % arithmetic, normalization, smoothing, etc., 
    % and plotting
    %
    % Petar Lambrev, 2012-2022

    properties (Dependent = true, SetAccess = private)
        dim % number of data points
    end
    
    methods
       % Constructor
       function SP = specdata(newX, newY, ID, varargin)
           % SPECDATA Create specdata objects
           %
           % Synthax
           % SP = specdata
           % SP = specdata(X, Y)
           % SP = specdata(X, Y, ID)
           % SP = specdata(X, Y, ID, ...)
           % SP = specdata(X, Y, ID, paramstruct)
           % SP = specdata(datastruct)
           % SP = specdata(specdata2Darray)
           %
           % Description
           % SP = specdata creates a blank spectrum SP
           % SP = specdata(X, Y) creates spectra with X and Y values
           % SP = specdata(X, Y, ID, ...) creates a spectrum SP with 
           %   X and Y values, ID, and additional parameters, 
           %   e.g. 'YType','CD','YUnit','mdeg',...
           % SP = specdata(datastruct) creates a spectrum or an array of
           %   spectra from datastruct
           % SP = specdata(specdata2Darray) converts specdata2D array to 
           %   specdata array
           %
           % Input parameters
           %   X - vector array
           %   Y - vector or matrix. If Y is a matrix, specdata(X, Y)
           %       creates an array of spectra for each column in Y
           %   ID - character array or cell array of characters
           %   paramstruct - struct of specdata parameters, e.g. 'YType'
           %   datastruct - struct or struct array with fields X,Y,ID,...
           %
           % Example
           %
           % x = (0:0.1:2*pi)';     % x: column vector
           % y = [sin(x), cos(x)];  % y: matrix with 2 columns
           % id = {'sinx', 'cosx'}; % id: array of 2 strings
           % A = specdata(x,y,id);  % A: array of 2 spectra
           
           proplist = {'DateTime','ExpID','XType','XUnit','YType','YUnit','Comment','History','Metadata'};           
           if exist('newX','var') || exist('newY','var')
               if istable(newX), newX = table2struct(newX); end
               if isa(newX, 'specdata2D')
                   % Convert specdata2D to specdata
                   nrows = numel(newX);
                   dim = newX(1).dim;
                   SP = specdata.empty(0,0);
                   for k = 1:nrows
                       if ~isequal(newX(k).dim, dim)
                           error('Cannot convert specdata2D object. The dimensions of the 2D spectra in the array are not equal.')
                       end
                       SP1 = specdata.empty(0,dim(2));
                       for col = 1:dim(2)
                           SP1(1,col).X = newX(k).X;
                           SP1(1,col).Y = newX(k).Y(:,col);
                           SP1(1,col).ID = [newX(k).ID, ' ', num2str(newX(k).T(col)), ' ', newX(k).TUnit];
                           for p = 1:numel(proplist)
                                   SP1(1,col).(proplist{p}) = newX(k).(proplist{p});
                           end
                       end
                       SP = [SP; SP1];
                   end
                   
               elseif isstruct(newX)
                   % Initialize from struct / table
                   n = numel(newX);
                   SP = specdata.empty(n,0);
                   for k = 1:n
                       newfields = fields(newX(k));
                       for field = newfields'
                           f = field{1};
                           if ismember(f, {'X', 'Y'})
                               if isvector(newX(k).(f)) && isnumeric(newX(k).(f))
                                   SP(k).(f) = newX(k).(f);
                               else
                                   error('X and Y must be numeric vectors.');
                               end
                           elseif strcmp(f, 'ID')
                               if isstring(newX(k).(f)) || (ischar(newX(k).(f)) && ~iscell(newX(k).(f)))
                                   SP(k).ID = newX(k).ID;
                               else
                                   error('ID must be a string.')
                               end
                           elseif ismember(f, proplist)
                               SP(k).(f) = newX(k).(f);
                           end
                       end
                   end
                   SP = reshape(SP, size(newX));
                   
               elseif isequal(size(newX),size(newY))
                   % Create from X, Y values
                   SP.X = newX(:);
                   SP.Y = newY(:);
                   if exist('ID','var')
                       SP.ID = ID;
                   end
                   numY = 1;
               else
                   X = reshape(newX,numel(newX),1);                   
                   if size(newY,1) == length(X)
                       numY = size(newY,2);
                   else
                       if size(newY,2) == length(X)
                          numY = size(newY,1);
                          newY = newY';  
                       else
                           error('the input X and Y arrays have incompatible sizes');
                       end
                   end
                   for j = 1:numY
                       if exist('ID','var')
                           if iscell(ID)
                               SP(j) = specdata(X,newY(:,j),ID{j});
                           else
                               if numY > 1
                                   SP(j) = specdata(X,newY(:,j),sprintf("%s %d",ID,j));
                               else
                                   SP(j) = specdata(X,newY(:,j),ID);
                               end
                           end
                       else
                           SP(j) = specdata(X,newY(:,j));
                       end
                   end
               end
           else
               SP.X = [];
               SP.Y = [];
           end
           if nargin > 3
               if isstruct(varargin(1))
                   optlist = varargin(1);
               else
                   optlist = struct(varargin{:});
               end
           else
               optlist = [];
           end
           
           if (exist('optlist','var') && isstruct(optlist))
               for j = 1:numY
                   for i = 1:length(proplist)
                       if(isfield(optlist,proplist{i}))
                           SP(j).(proplist{i}) = optlist.(proplist{i});
                       end
                   end
               end %for
           end %if
        end %Constructor
               
       %% Queries
       function n = get.dim(SP)
           % dim - number of data points
           %
           % Synthax
           %    n = MySpectrum.dim
           n = zeros(1,length(SP));
           for i = 1:length(SP)
               n(i) = length(SP(i).X);
           end %for
       end %size

       function res = mldivide(A,B)
           % MLDIVIDE Solve systems of linear equations Ax = B for x
           % 
           % Synthax
           % x = A \ B
           % x = mldivide(A,B)
           %
           % Description
           % x = A\B solves the system of linear equations A*x = B,
           % where A and B are spectra and x is a scalar
           %
           % Examples
           %
           % Estimate concentration using standard spectrum.
           % Let A is the normalized reference absorption spectrum of a
           % compound and B is a measured spectrum of a solution of the
           % same compound. The concentration of the solution is:
           % 
           % c = A \ B;
           %
           % If the result is correct, the measured spectrum should be 
           % identical to the reference multiplied by concentration:
           % 
           % plot(A*c, B);
           %
           % Find the composition of a mixture.
           % Let A is an array of reference spectra for N different compounds.
           % Let B is a measured spectrum of a mixture of the N compounds.
           % To find the composition of the mixture, execute:
           % 
           % c = A \ B;
           %
           % The result c is an N-element array satisfying the equation:
           %
           % B = A(1)*c(1) + A(2)*c(2) + ...
           %
           % See also MLDIVIDE
           
           Amat = A.xymat; Amat(:,1) = [];
           Bmat = B.xymat; Bmat(:,1) = [];
           res = Amat\Bmat;               
       end %rdivide
              
       function [TF,P] = islocalmax(SP,varargin)
           % ISLOCALMAX Find local maxima
           %
           % Synthax
           %   [TF,P] = islocalmax(SP)
           %   [TF,P] = islocalmax(SP,Name,Value)
           %
           % Applies the namesake MATLAB function to the spectra
           %
           % See also: ISLOCALMIN

           [TF,P] = islocalmax([SP(:).Y],varargin{:});
       end

       function [TF,P] = islocalmin(SP,varargin)
           % ISLOCALMIN Find local minima
           %
           % Synthax
           %   [TF,P] = islocalmin(SP)
           %   [TF,P] = islocalmin(SP,Name,Value)
           %
           % Applies the namesake MATLAB function to the spectra
           %
           % See also: ISLOCALMIN

           [TF,P] = islocalmin([SP(:).Y],varargin{:});
       end
       
        function res = xymat(SP)
           % XYMAT Matrix of X-Y values
           %
           % Synthax
           %    res = Sp.xymat;
           %
           % Description
           % The function returns a matrix res with m rows and n+1 columns,
           % where m is the number of X values and n is the number of spectra in Sp.
           % The first column contains the X data.
           
           l = length(SP);
           res = zeros(SP(1).dim,l+1);
           res(:,1) = SP(1).X;
           for i = 1:l
               if ~isequal(SP(i).X, SP(1).X)
                   error('Cannot create XY matrix. X data do not match.');
               end
               res(:,i+1) = SP(i).Y;
           end           
        end      
       
        function tbl = xytable(SP)
            % XYTABLE Create an X-Y table from specdata objects
            %
            % Synthax
            %    tbl = xytable(SP)
            %
            % Description
            % The spectra in SP must have a common X axis.
            % The output tbl is a table with the X axis in the first column
            % and Y values for all spectra in subsequent columns with the
            % corresponding IDs as variable names (column headings).   
            % 
            % The IDs are modified to valid MATLAB identifiers
            % (spaces and special symbols are removed).
            %
            % The original IDs are saved as variable descriptions
            % (accessed as tbl.Properties.VariableDescriptions).
            
            n = numel(SP);
            a = zeros(SP(1).dim, n+1);
            a(:,1:2) = [SP(1).X, SP(1).Y];
            if n > 1
                for k = 2:n
                    if comp(SP(1),SP(k))
                        error('All spectra must have the same X axis. Use setX to make the X axes equal.');
                    end
                    a(:,k+1) = SP(k).Y;
                end
            end
            vardescr = [{SP(1).XType}, {SP.ID}];            
            varnames = matlab.lang.makeValidName(vardescr, 'ReplacementStyle','underscore');
            varnames = matlab.lang.makeUniqueStrings(varnames);
            tbl = array2table(a, 'VariableNames', varnames);            
            tbl.Properties.VariableUnits = [{SP(1).XUnit}, {SP.YUnit}];
            tbl.Properties.VariableDescriptions = vardescr;
        end
           
       %% Other calculations
       function res = baseline(SP,X1,X2)
           % Subtract sloped baseline based on two points
           %
           % Syntax
           %  NewSpectrum = MySpectrum.baseline(X1,X2)
           %
           % Description
           % Subtracts straight baseline that is calculated as a line
           % crossing two wavelengths (X1, X2) in the original spectra.
           % The resulting spectra will have zero amplitude at X1 and X2.
           
           res = SP;
           for i = 1:length(SP)
               P1 = SP(i).Yx(X1);
               P2 = SP(i).Yx(X2);
               baseX = SP(i).X;
               YC = (P2-P1)/(X2-X1);
               baseY = P1 + (baseX-X1).*YC;
               res(i).Y = SP(i).Y - baseY;
           end
       end
       
       function [fitres, gof] = fit(SP,model,varargin)
            % FIT Fit model to data
            %
            % Syntax
            %     fitres = SP.fit(model, varargin)
            %
            % Description 
            % Fits the spectrum SP using the specified model.
            %
            % The model can be either a built-in library model, e.g.
            % 'poly5' - 5th-order polynomial or
            % 'gauss5' - sum of 5 gaussians, or
            % it can be a custom nonlinear model defined using fittype.
            %
            % Specifics using a gaussian model ('gauss1', 'gauss2', ...):
            %    By default the following constraints are applied:
            %    Positive amplitudes;
            %    Positive center X values (wavelengths);
            %    Center X values fall within the X data limits;
            %    Widths cannot be larger than the total X data range.
            %
            % The returned result fitres is an object of class cfit.
            % A second optional returned result gof is the goodness of fit.
            %
            % The function requires MATLAB Curve-Fitting Toolbox.
            %
            % Note: Due to a limitation in MATLAB's fit function, only
            % one spectrum is fitted. Additional data in SP are ignored.
            % Subject to change in a future versioni.
            %
            % See also:
            %    fit, cfit, fittype
            
            
            % Set model type
            if length(SP)>1
                SP = SP(1);
                warning('Warning. Only first spectrum will be fitted.');
            end
            if ischar(model)
                fmodel = fittype(model);
                modelname = model;
            else
                if isa(model,'fittype')
                    fmodel = model;
                    modelname = type(model);
                else 
                    error 'The argument model must be either a string specifying library model or a fittype object.';
                end
            end
            if nargin == 3
                options = varargin{1};
            else
                options = fitoptions(fmodel);
            end
            if nargin > 3
                set(options,varargin{:});
            end
            
            % gauss model
            if contains(modelname,'gauss') && nargin ~= 3
                ncomp = str2num(modelname(6)); % number of components
                defopt = fitoptions(fmodel);
                
                dat = SP(1);
                % set lower bounds
                LB = get(options,'Lower');
                if isequal(LB,get(defopt,'Lower'))
                    minX = min(SP.X);
                    for i = 0:ncomp-1
                        LB(i*3+1) = 0; % positive amplitude
                        LB(i*3+2) = minX; % minimum center
                    end
                    set(options,'Lower',LB);
                end
                UB = get(options,'Upper');
                if isequal(LB,get(defopt,'Lower'))
                    minX = min(SP.X);
                    maxX = max(SP.X);
                    rngX = maxX-minX;
                    maxY = max(SP.Y);
                    for i = 0:ncomp-1
                        UB(i*3+1) = maxY; % maximum amplitude
                        UB(i*3+2) = maxX; % maximum center
                        UB(i*3+3) = rngX; % maximum width
                    end
                    set(options,'Upper',UB);
                end
            end          
            [fitres, gof] = fit(SP.X, SP.Y, model, options);
        end %fit
        
       function res = merge(SP, weights)
            % MERGE Join spectra into one
            % 
            % Synthax
            %     B  = A.merge
            %     B  = A.merge(weights)
            %     B  = merge(A,...)            
            %
            % Description
            % All spectra in A will be merged into a single spectrum B
            % The resulting object B contains all X values found in SP
            % Y values in B that correspond to the same X value are averaged
            %
            % Note: the function is similar to mean and returns the same
            % result if all spectra in A share the same X values.
            % However if the X values are different, mean returns an error,
            % whereas merge returns a combined spectrum with averaged
            % overlaps. Therefore merge is useful for combining
            % measurements done in different spectral (wavelength, etc.)
            % ranges. 
            %
            % B = merge(A,weights) joins by using weighted average
            %
            % Input parameters
            %   A - array of spectra (specdata)
            %   weights (vector) statistical weight for each spectrum            
            
            res = SP(1);
            if length(SP) <=1
                 return;
            end
            
            if ~exist('weights','var')
                weights = ones(length(SP));
            end
            if length(weights) < length(SP)
                error('The number of weight factors must equal the number of spectra to be merged.');
            end
            
            % combine all X values into one array
            Xall = SP(1).X;
            for i=2:length(SP)
                Xall = union(Xall,SP(i).X);
            end
            
            % get the unique X values
            X = unique(Xall);
            res.X = [];
            res.X = X;
            res.Y = zeros(length(X),1);
            
            % for each unique X value, average all Y values
            for j = 1:length(X)
               Ynum = 0; Ysum = 0;
               for i = 1:length(SP)
                   xi = find(SP(i).X==X(j),1);
                   if ~isempty(xi)
                       Ysum = Ysum + weights(i)*SP(i).Y(xi);
                       Ynum = Ynum + weights(i);
                   end %i
               end %j
               res.Y(j,1) = Ysum / Ynum;
            end
        end
            
       %% Display and plotting
              function p = plot(SP,options)
           % PLOT Plot spectra
           %
           % Synthax
           %    plot(A)
           %    p = plot(A,Name,Value)
           %
           % Description
           % plot(A) plots the array of spectra A in the current figure. 
           % If there is no figure, creates a new one. 
           % The axes are labeled using XType/XUnit, YType/YUnit of the 
           % first spectrum in A. Legend shows spectra IDs.
           %           
           % Name-Value arguments
           %    LegendText - list of properties to display as legend text 
           %         (default is "ID", can also be custom metadata)
           %    LegendFun - custom legend function. Example: 
           %         plot(A, "LegendFun", @(x) x.ExpID+"/"+x.ID)
           %    LegendInterpreter - legend interpreter ("tex","none")
           %    LegendBox - legend box ("on","off")
           %    LegendLocation - legend location ("best")
           %    SmoothLine - spline interpolation between data points
           %    (true,false)
           %    Line properties (LineWidth, Marker, etc.)
           
           arguments
               SP specdata {mustBeNonempty}
               options.LegendText string = "ID"
               options.LegendFun function_handle = function_handle.empty;
               options.LegendInterpreter string = "tex"
               options.LegendLocation string = "best"
               options.LegendBox string = "off"
               options.SmoothLine (1,1) logical = false
               options.?matlab.graphics.chart.primitive.Line
           end
           
           SP = SP(:);
           plotOptions = rmfield(options,["LegendText","LegendFun","LegendBox","LegendInterpreter","LegendLocation","SmoothLine"]);
           Xmin = min(SP(1).X); Xmax = max(SP(1).X);
           
           h = matlab.graphics.chart.primitive.Line.empty;
           for i = 1:length(SP)
               % set Xmin and Xmax for Xlim
               Xmin_i = min(SP(i).X);
               Xmax_i = max(SP(i).X);
               if Xmin_i < Xmin
                   Xmin = Xmin_i;
               end
               if Xmax_i > Xmax
                   Xmax = Xmax_i;
               end
               % plot curve
               nx = numel(SP(i).X) ;
               if options.SmoothLine
                   nx = (nx-1) * 10;
                   dx = (max(SP(i).X) - min(SP(i).X)) / nx;
                   try
                       smx = min(SP(i).X):dx:max(SP(i).X);
                       smy = interp1(SP(i).X,SP(i).Y,smx,'spline');
                       handle = plot(smx,smy,'-o','MarkerIndices',1:10:numel(smx),plotOptions);
                   catch
                       warning('Cannot interpolate data for smoothing.')
                       handle = plot(SP(i).X,SP(i).Y,plotOptions);                       
                   end
               else
                   handle = plot(SP(i).X,SP(i).Y,plotOptions);
                   
               end
               h = [h handle];
               hold all;
           end
           
           % Create Legend
           if ~isempty(options.LegendFun)
               LString = options.LegendFun(SP);
           else
               if isempty(options.LegendText)
                   options.LegendText = "ID";
               end
               LText = options.LegendText(:);
               PropTable = SP.proptable;
               if isempty(setdiff(LText,PropTable.Properties.VariableNames))
                   % Display specdata properties
                   LString = join(string(PropTable{:,LText}),2);
               else
                   % Display custom legend text
                   LString = LText;
               end
           end
           lh = legend(h,LString,"Location",options.LegendLocation);
           lh.Interpreter = options.LegendInterpreter;
           lh.Box = options.LegendBox;
           
           if strcmp(SP(1).XUnit,'')
               Lbl = sprintf('%s',SP(1).XType);
           else
               Lbl = sprintf('%s (%s)',SP(1).XType,SP(1).XUnit);
           end
           xlabel(Lbl);
           if strcmp(SP(1).YUnit,'')
               Lbl = sprintf('%s',SP(1).YType);
           else
               Lbl = sprintf('%s (%s)',SP(1).YType,SP(1).YUnit);
           end
           xlim([Xmin Xmax]);
           ylabel(Lbl);

           
           if nargout, p = h; end
       end %plot
       
       function ploterror(S,E,varargin)
           % PLOTERROR Plot spectra with shaded error
           %
           % Synthax
           %   plot(S,E)
           %   plot(S,E,Name,Value)
           %
           % Description
           %   plot(S,E) plots the spectra S as lines and the associated
           %   error intervals E as a shaded area around them. S and E are
           %   specdata objects

           % Get color order
           col = get(gca,'ColorOrder');
           hold on
    
           % Plot the shaded area first
           for k = 1:numel(S)
               xx = [S(k).X; S(k).X(end:-1:1)];
               yy = [S(k).Y + E(k).Y; S(k).Y(end:-1:1) - E(k).Y(end:-1:1)];
               p = fill(xx,yy,col(k,:),'HandleVisibility','off');
               p.EdgeColor = 'none';     alpha(p,0.25)
           end

           % Plot the central lines
           set(gca,'ColorOrderIndex',1)

           if isempty(varargin)
                plot(S);
           else
               plot(S,varargin{:})
           end
       end

       function w = fwhm(S)
           %FWHM Find Full Width at Half Maximum of spectra
           %   Synthax
           %     w = fwhm(S)
           %
           %   Description
           %     w = fwhm(S) returns the FWHM of the spectra S. Uses spline interpolation
           %     to fill in the missing points.
           
           n = numel(S); % number of spectra
           w = zeros(n,1);
           
           for k = 1:n
               % Spline interpolation
               
               % Create fine grid
               xi = linspace(min(S(k).X),max(S(k).X),numel(S(k).X)*100)';
               
               % interpolate
               yi = interp1(S(k).X,S(k).Y,xi,'spline');
               
               % find maximum
               [ymax, imax] = max(yi);
               hmax = ymax/2;
               
               % find left and right half-max
               [~, ileft] = min(abs(hmax-yi(1:imax)));
               [~, iright] = min(abs(hmax-yi(imax:end)));
               
               % fwhm is the difference
               w(k) = abs(xi(ileft)-xi(imax+iright-1));
           end % fwhm
       end
           
       function markpeaks(spect,delta,varargin)
           % MARKPEAKS Mark peaks on the current graph
           % 
           % Synthax
           %    markpeaks(A)  
           %    markpeaks(A, delta)
           %    markpeaks(A, delta, [Text properties])
           %
           % Description
           % markpeaks(A) detects peaks (local minima and maxima) in the 
           % spectrum A and writes their position (X) on the current
           % figure.
           %
           % markpeaks(A, delta) marks only peaks with prominence higher
           % than delta.
           %
           % Input parameters
           % A - specdata or specdata array
           % delta - minimal amplitude of the peaks (along the Y axis)
           %      to be marked. If not supplied, delta is assumed to be half
           %      the range of Y axis values.
           % [Text parameters] - optional arguments specifying text
           %      properties, e.g. 'Color', 'FontSize', etc.
           %
           % Example
           %
           % Chl = [Chla, Chlb];
           % delta = max(Chl) / 2;
           % figure; plot(Chl) % plot two spectra
           % markpeaks(Chl, delta) % mark detected peaks
           %
           % figure; plot(Chl)
           % markpeaks(Chl, [436 663]) % mark selected peaks for both
           %
           % figure; plot(Chl)
           % markpeaks(Chl, {[436 663],[480 652]}) % mark different peaks
           %
           % see also: FINDLOCALMAX, FINDLOCALMIN
           
           n = length(spect);
           co = get(gca,'ColorOrder'); % color scheme
           conum = size(co,1); % number of colors
           hold on
           for spi = 1:n
               Y = spect(spi).Y;
               X = spect(spi).X;
               if (nargin < 2)
                   delta = (max(Y) - min(Y)) / 4;
               end

               % Detect peak positions
               tmax = find(islocalmax(Y,'MinProminence',delta));
               tmin = find(islocalmin(Y,'MinProminence',delta));
               l = ylim; % xl = xlim;
               d = (l(2)-l(1))*0.04;

               % get current color index
               coi = rem(spi,conum); if ~coi, coi = conum; end
               if any(tmax)
                   for i = 1:height(tmax)
                       px = X(tmax(i)); py = Y(tmax(i)) + d;
                       if (py<l(2))
                           text(px,py,sprintf('%.0f',X(tmax(i))),'HorizontalAlignment','center','color',co(coi,:), varargin{:});
                       end
                   end
               end
               if any(tmin)
                   for i = 1:height(tmin)
                       px = X(tmin(i)); py = Y(tmin(i)) - d;
                       if py>l(1)
                           text(px,py,sprintf('%.0f',X(tmin(i))),'HorizontalAlignment','center','color',co(coi,:),varargin{:});
                       end
                   end
               end
           end
       end %markpeaks
       
       %% Save
       function save(SP,filename)
           % SAVE Save spectra as tab-delimited ASCII file.
           %
           % Synthax
           %    A.save(filename)
           %
           % Description
           % If A is an array of spectra, they are saved in the same file 
           % as separate columns. 
           % 
           % To save the spectra as a MATLAB file (.mat) instead, 
           % use the MATLAB built-in save command:
           %    save filename A 
           
           outheader = strcat(SP(1).XType,'\t',SP(1).ID);
           outnum = cat(2,SP(1).X,SP(1).Y);
           if length(SP) > 1               
               for i = 2:length(SP)
                   if ~isequal(SP(i).X, SP(1).X)
                       error('To save multiple spectra in a single ASCII file, their X axes must be identical. If needed, use the setX method to convert the axes.');
                   end
                   outheader = strcat(outheader,'\t',SP(i).ID);
                   outnum = cat(2,outnum,SP(i).Y);
               end
           end
           outheader = strcat(outheader,'\r\n');
           fid = fopen(filename,'w');
           fprintf(fid,outheader);
           for i = 1:size(outnum,1)
               ln = sprintf('%g',outnum(i,1));
               for j = 2:size(outnum,2)
                   ln = sprintf('%s\t%g',ln,outnum(i,j));
               end
               ln = strcat(ln,'\r\n');
               fprintf(fid,ln);
           end
           fclose(fid);
          fprintf('File %s saved.\r',filename);           
       end %save
       
       function saveh5(SP,FileName)
           % SAVEH5 Save spectra as an HDF5 file.
           %
           % Synthax
           %    A.saveh5(filename)
           
           Profile = 'Simple';
           % Get properties
           switch Profile
           case 'Simple'
                   Props = {'XType', 'XUnit', 'YType', 'YUnit'};
               otherwise
                   Exclude = {'X','Y','ID'};
                   AllProps = properties(SP(1));
                   Props = AllProps(~ismember(AllProps,Exclude));
           end               
           
           if exist(FileName,'file')
               delete(FileName)
           end
           
           for k = 1:numel(SP)
              DS = ['/' SP(k).ID];              
              h5create(FileName, DS, [2 SP(k).dim])
              h5write(FileName, DS, [SP(k).X'; SP(k).Y'])
              for p = 1:length(Props)                  
                  h5writeatt(FileName, DS, Props{p}, SP(k).(Props{p}))
              end
           end
           
           h5writeatt(FileName,'/','Creator','spectromatic v2.3');
           h5writeatt(FileName,'/','Profile',Profile)
           
           fprintf('File %s saved.\r',FileName);
       end %saveh5       

       function data2D = unstack(data, expr)
           % UNSTACK Convert specdata array to specdata2D array
           %
           % Synthax
           %   data2D = unstack(data)
           %   data2D = unstack(data, expr)
           %
           % Input/Output parameters
           %   data - specdata array
           %   expr - regular expression that defines the T dimension
           %   data2D - specdata2D array
           %                       
           
           id1 = get(data,'ID');
           id2 = regexprep(id1, expr, "");
           id3 = regexp(id1, expr, 'match', 'once');
           idn = str2double(regexp(id3, '\d+', 'match', 'once'));
           
           [g, gid] = findgroups(id2);
           n1 = numel(gid);
           
%            if n1 < 2
%                error('The data could not be divided into groups.');
%            end
           
           data2D = specdata2D.empty(n1,0);
           
           for k = 1:n1
               l = find(g==k);
               data2D(k,1) = data(l)';
               data2D(k,1).ID = gid{k};
               data2D(k,1).T = idn(l)';
           end           
       end % unstack
       
    end %methods
    
    %% Private methods
    methods (Access = protected, Hidden = true)
        function match = comp(SP1, SP2)
            % compare two spectra if their X dimensions match
            % return 1 if there is a difference and 0 if there is no diff
            if isequal(SP1.X,SP2.X)
                match = 0;
            else
                match = 1;
            end
        end
    end % private methods
    
    %% Load data
    methods (Static)
        function S = load(FilePattern,varargin)
            % LOAD Load spectra from ASCII file
            %
            % Synthax
            %    data = specdata.load(files,[options])
            %
            %    files: a cell array of filenames, 
            %       e.g. {'file1.txt';'file2.txt';'file3.txt'}
            %       or a string containing wildcards, e.g. '*.txt'
            %
            %    Optional parameters:
            %    - FileType - data file format. Can be either:
            %           'XY' - tab- or space-delimited numeric data      
            %           'ASCII' - generic text file, read numeric data only                      
            %           'JWS' - Jasco text file
            %           'Chirascan' - ChiraScan (ProData) CSV or TXT file
            %           'SpectraSuite' - Ocean Optics SpectraSuite
            %           'MultiCol' - multiple Y columns and header
            %           (each column will be treated as a spectrum)
            %           'spreadsheet' - CSV or Excel spreadsheet
            %
            %    If FileType is not specified, MATLAB version 2016b will
            %    load data from a text file, CSV file or a spreadsheet, 
            %    using the readtable function. 
            %    In earlier MATLAB versions only numeric data will be
            %    loaded. Lines containing other characters will be ignored.
            %
            %    - Delimiter - specify column delimiter for auto and
            %            MultiCol file formats
            %
            %    - Decimal - specify decimal separator, e.g. ','. Only applies
            %      to ASCII file type
            %
            %    - XType, XUnit, YType, YUnit, ExpID, Comment
            %
            %    Optional parameters can be passed as separate arguments, e.g.            %
            %       data = specdata.load('file1.txt','Type','Abs','XUnit','nm')
            %
            %    Alternatively, all optional parameters can be given
            %    in a single structure variable:            %
            %       options.XType = 'Wavelength';
            %       options.XUnit = 'nm';
            %       options.FileType = 'XY';
            %       data = specdata.load('file1.txt',options);
            FilePattern = string(FilePattern);
            if numel(FilePattern)>1
                Files = FilePattern;
            else
                Files = getdir(FilePattern, 'FullPath', true);
            end %if
            if isempty(Files) || length(Files)<1, error('File not found.'); end
            disp('Loading spectra...');
            if nargin > 1
                if isstruct(varargin(1))
                    Options = varargin(1);
                else
                    Options = struct(varargin{:});
                end
            else
                Options = [];
            end                
            S = specdata.empty(length(Files),0);
            multicol = 0;
            n = 0;
            for j = 1:length(Files)
                if isfield(Options,'Delimiter')
                    delim = Options.Delimiter;
                else                    
                    delim = -1;
                end
                if isfield(Options,'Variable')
                    % Variable (YType) to load from data file
                    % Applicable to Chirascan
                    loadvariable = Options.Variable;
                else
                    loadvariable = '';
                end
                if isfield(Options,'FileType')
                    switch Options.FileType 
                        case 'XY' 
                            d = specdata.loadXY(Files{j});
                        case 'JWS'
                            [d, DateTime, YType, Comment] = specdata.loadJWS(Files{j}); 
                            XType = 'Wavelength'; xunit = 'nm';
                        case {'ASCII','ascii'}
                            if isfield(Options,'Decimal')
                                DecimalSeparator = Options.Decimal;
                                d = specdata.loadascii(Files{j},delim,DecimalSeparator);
                            else
                                d = specdata.loadascii(Files{j},delim);
                            end
                        case {'Chirascan','chirascan'}
                            [d, DateTime,XType, YType, Comment] = specdata.loadChirascan(Files{j},loadvariable);
                        case 'SpectraSuite'
                            [d, DateTime] = specdata.loadSpectraSuite(Files{j});                              
                        case {'MultiCol','Multicol','MultiColumn','Multicolumn'}
                            multicol = 1;
                            [d, Comment] = specdata.loadMultiCol(Files{j},delim);
                        case {'SpreadSheet','spreadsheet','text'}
                            [d, XType, YType] = specdata.loadSpreadsheet(Files{j},varargin);
                        otherwise                            
                            error('%s is not a recognized file type.',Options.FileType);
                    end
                else
                    % Get MATLAB Version
                    MATLAB_Version = ver('MATLAB');
                    MATLAB_Version = str2double(MATLAB_Version.Version);
                    if MATLAB_Version >= 9.1   
                        [d, XType, YType] = specdata.loadSpreadsheet(Files{j},varargin);
                    else
                        d = specdata.loadascii(Files{j},delim);
                    end
                end
                d = sortrows(d,1);
                
                for k = 1:(size(d,2))-1
                    S(j,k) = specdata(d(:,1),d(:,k+1));
                    if exist('Comment','var')
                        S(j,k).Comment = string(Comment);
                    end
                    S(j,k).History = sprintf("load %s",Files{j});
                    if exist('YType','var')
                        if iscell(YType)                            
                            YUnit = regexp(YType{k},'\w*','match');
                            if ~isempty(YUnit)
                                S(j,k).YType = string(YUnit{1});
                                if length(YUnit)>1
                                    S(j,k).YUnit = string(YUnit{2});
                                end
                            end
                        else
                            S(j,k).YType = string(YType);
                        end
                    end
                    if exist('XType','var')
                        if iscell(XType)
                            S(j,k).XType = string(XType{k});
                        else
                            S(j,k).XType = string(XType);
                        end
                    end
                    if exist('xunit','var')
                        if iscell(XType)
                            S(j,k).XUnit = string(xunit{k});
                        else
                            S(j,k).XUnit = string(xunit);
                        end
                    end                    
                    ft = dir(Files{j});
                    if (~exist('datetime','var'))
                        DateTime = ft.datenum;
                    end
                    if isdatetime(DateTime)
                        S(j,k).DateTime = DateTime;
                    elseif isnumeric(DateTime)
                        S(j,k).DateTime = datetime(DateTime,'ConvertFrom','datenum');
                    else
                        S(j,k).DateTime = datetime(DateTime,'InputFormat','yyyy-MM-dd HH:mm:ss');
                    end
                    [DirName, FileName] = fileparts(Files{j});
                    % get relative path
                    s = strfind(DirName, filesep);
                    if ~isempty(s), DirName = DirName(s(end)+1:end); end                    
                    S(j,k).ExpID = string(DirName);
                    S(j,k).ID = string(regexprep(FileName, "\.(\w{3})$",""));
                    
                    props = {'XType','XUnit','YType','YUnit','ExpID','Comment','Date'};
                    if (exist('Options','var'))
                        for p = 1:length(props)
                            if (isfield(Options,props{p}))
                                prop = Options.(props{p});
                                if iscell(prop)
                                    S(j,k).(props{p}) = string(prop{j});
                                else
                                    S(j,k).(props{p}) = string(prop);
                                end
                            end
                        end %for p
                    end %if exist
                end %for k
                fprintf('File %d: %s\n',j,Files{j});
            end %for j
        end       
        
    end %Static
    
    methods (Static, Access = private, Hidden = true)
        function d = loadXY(filename)
            d = load('-ascii',filename);
            if (size(d,2) < 2), error('File must contain two columns of numerical data'); end
        end %loadXY
        
        function d = loadascii(FileName,Delimiter,DecimalSeparator)
            fi = fopen(FileName);
            if (fi == -1)
                error('File %s not found.', FileName);
            end %if
            [~,~,ext] = fileparts(FileName);
            if strcmp(ext,'.csv') && Delimiter == -1
                Delimiter = ',';
            end
            
            j = 0;
            d = zeros(2000,2);
            while ~feof(fi)
                l = fgetl(fi);                
                if isempty(l) || ~isstrprop(l(1),'digit')
                    if j == 0, continue
                    else, break
                    end
                end
                if exist('DecimalSeparator','var')
                    l = strrep(l,DecimalSeparator,'.');
                end
                if Delimiter ~= -1
                    t = textscan(l,'%f','delimiter',Delimiter);
                else
                    t = textscan(l,'%f');
                end
                nt = length(t{1}) ;
                if nt < 2, continue, end
                if nt > size(d,2)
                    d(:,nt) = zeros(size(d,1),1);
                end
                j = j + 1;                
                d(j,:) = t{1};
            end %while
            fclose(fi);
            d(j+1:end,:) = [];
            if (size(d,1) < 1), error('File could not be read. No numerical data found.'); end
        end %loadascii
        
        function [d,header] = loadMultiCol(f,delim)
            fi = fopen(f);
            if (fi == -1)
                error('File %s not found.', f);
            end %if
            j = 0;
            d = zeros(1000,2);
            numcol = 0;
            header = {};
            while ~feof(fi)
                l = fgetl(fi);                
                if isempty(l), continue, end
                if isstrprop(l(1),'digit')
                    %numeric data -> read data;
                    if numcol == 0
                        if delim ~= -1
                            t = textscan(l,'%f','delimiter',delim);
                        else
                            t = textscan(l,'%f');
                        end
                        numcol = length(t{1});
                        d = zeros(1000,numcol);
                    else
                        if delim ~= -1
                            t = textscan(l,'%f',numcol,'delimiter',delim);
                        else
                            t = textscan(l,'%f',numcol);
                        end
                    end
                    if length(t{1}) < 2, continue, end
                    j = j + 1;
                    d(j,:) = t{1}; 
                else
                    %non-numeric data -> read header;
                    t = textscan(l,'%s');
                    header = t{1};
                    numcol = length(header);
                    d = zeros(1000,numcol);
                end %if                
            end %while
            
            fclose(fi);
            d(j+1:end,:) = [];
            if (size(d,1) < 1), error('File could not be read. No numerical data found.'); end
        end %loadMultiCol

        function [d, datetime, xunits, yunits, comment] = loadChirascan(filename,loadvariable)
            fi = fopen(filename);
            if (fi == -1)
                error('File %s could not be opened.', filename);
            end %if
            xunits = {}; yunits = {};
            comment = '';
            datetime = '';
            begindata = 0;
            ncol = 1;
            % Read Header
            l = fgetl(fi);
            if ~strcmp(l(1:7),'ProData')
                error('File does not appear to have a valid Chirascan format.');
            end
            while begindata==0
                l = fgetl(fi);
                if strcmpi(l,'Data:')
                    l = fgetl(fi);
                    begindata = 1;
                else
                    comment = [comment, newline, l];
                    if strfind(l,'Available Properties:')
                        s = regexp(l,'\d+$','match');
                        ncol = str2double(s{1});
                        for pl = 1:ncol
                            l2 = fgetl(fi);
                            s = regexp(l2,'^\w*\>','match');
                            props(pl) = s(1);
                        end
                    elseif strfind(l,'Available Dimensions:')
                        s = regexp(l,'\d+$','match');
                        ndim = str2double(s{1});
                        if ndim > 2
                            error('Unsupported Chirascan file format.');
                        end                        
                        l2 = fgetl(fi);                        
                        xunit = regexp(l2, '^\w*\>','match'); xunit = xunit{1};                        
                    end                    
                end
            end
            % Read Data
            if ~isempty(loadvariable)
                [varfound, varindex] = ismember(loadvariable, props);
                if ~varfound
                    varindex = []; d = [0 0]; xunits = {''}; yunits = {''};
                    warning('Variable %s not found in data file',loadvariable)
                    fclose(fi);
                    return
                end
            else
                varindex = 1:ncol;
            end

            if ndim==1
                % read one dimension data
                for k = 1:ncol
                    xunits{k} = ''; yunits{k} = '';
                    l2 = fgetl(fi); l3 = fgetl(fi);
                    s = textscan(fi,'%f,%f');
                    if isempty(loadvariable)                        
                        xunits{k} = l2; yunits{k} = l3;
                        d(:,1) = s{1}; d(:,k+1) = s{2};
                    else
                        if strcmp(l3, props{varindex})
                            xunits = l2; yunits = l3;
                            d(:,1) = s{1}; d(:,2) = s{2};
                            break
                        end
                    end
                    if feof(fi), break; end
                end
            else % ndim == 2
                % read two-dimension data
                d = [];
                for k = varindex
                    looking = true;
                    while looking || feof(fi)
                        l = fgetl(fi);
                        if strcmp(l,props{k})
                            yunits{k} = props{k};
                            xunits{k} = xunit;
                            looking = false;
                            l2 = fgetl(fi);
                            l3 = fgetl(fi);
                            firstdigit = true;
                            rownum = 1;
                            while firstdigit
                                l = fgetl(fi);
                                if regexp(l,'^\d')==1
                                    s = textscan(l,'%f','Delimiter',', \t');
                                    s = s{1};
                                    if isnumeric(s)                                        
                                        if numel(varindex)>1
                                            d(rownum,k+1) = mean(s(2:end));
                                            d(rownum,1) = s(1);
                                        else
                                            d(rownum,:) = s';
                                        end
                                        rownum = rownum + 1;
                                    end
                                else
                                    firstdigit = false;
                                end
                            end
                        end
                    end
                end
            end
            while ~feof(fi)
                l = fgetl(fi);
                if regexp(l,'^Time Stamp :')
                    datetime = strtrim(l(13:end));
                end
            end
            fclose(fi);
            if numel(varindex) == 1 && ndim == 2
                yunits = repmat(yunits(varindex),1,size(d,2));
                xunits = repmat(xunits(varindex),1,size(d,2));                
            end
            xunits = regexprep(xunits,',$','');
        end
        
        function [d, DateTime, YUnits, Comment] = loadJWS(filename)
            fi = fopen(filename);
            if (fi == -1)
                error('File %s could not be opened.', filename);
            end %if
            begindata = 0;
            YUnits = cell(0,3);
            Comment = '';
            while (begindata==0)
                l = fgetl(fi);
                if (l == -1), error('Invalid file format.'), end
                if (strfind(l,'YUNITS'))
                    ChannelNum = 1;
                    YUnits{1} = l(8:end);
                end                    
                if(strcmp(l(1:5),'TITLE'))
                    Comment = l(7:end);
                end                
                if (strfind(l,'Y2UNITS'))
                    ChannelNum = 2; 
                    YUnits{2} = l(8:end);
                end                
                if (strfind(l,'Y3UNITS'))
                    ChannelNum = 3; 
                    YUnits{3} = l(8:end);
                end
                if(strcmpi(l(1:4),'DATE'))      
                    try
                    datev = datevec(l(6:end),'yy/mm/dd');
                    catch
                        warning('Could not interpret Date in file %s. Using file date/time instead',filename)
                    end
                end
                if(strcmpi(l(1:4),'TIME'))
                    try
                    timev = datevec(l(6:end),'HH:MM:SS');
                    catch
                         warning('Could not interpret Time in file %s. Using file date/time instead',filename)
                    end
                end
                if (strcmpi(l,'XYDATA')), begindata = 1; end
            end
            switch ChannelNum
                case 1
                    ff = '%f %f';
                case 2
                    ff = '%f %f %f';
                case 3
                    ff = '%f %f %f %f';
            end
            d = textscan(fi,ff);
            Comment = [Comment,sprintf('\n')];
            while ~feof(fi)
                Comment = [Comment,fgets(fi)];
            end
            fclose(fi);
            d = cat(2,d{:});
            if (size(d,2) < 2), error('File must contain two columns of numerical data'); end
            % d = [d{:,1},d{:,2}];            
            
            if (exist('datev','var') && exist('timev','var'))                
                DateTime = datetime(cat(2,datev(1:3),timev(4:6)));
            else
                fd = dir(filename);
                DateTime = datestring(fd.datenum,'ConvertFrom','datenum');
            end
        end %loadJWS
        
        function [d, datetime] = loadSpectraSuite(filename)
            % load SpectraSuite data file
            fi = fopen(filename);
            if (fi == -1)
                error('File %s not found.', filename);
            end %if
            l = fgetl(fi); % first line
            if ~strfind(l,'SpectraSuite')
                error('%s is not a SpectraSuite data file.',filename);
            end
            begin = 0;
            
            % parse header
            while ~feof(fi)
                l = fgetl(fi); % second line
                if strcmp(l(1:5),'Date:')
                    datetime = l(7:end);
                    continue
                else
                    if strfind(l,'Begin')
                        begin = 1;
                        break
                    end
                end
            end
            
            if begin == 0
                error('No data found in %s',filename);
            end
            
            % parse data
            d = [];
            while ~feof(fi)
                l = fgetl(fi);
                if isempty(l) || ~isempty(strfind('l','End'))
                    break
                else
                    ls = textscan(l,'%f',2);
                    if ~isempty(ls{1})
                        d(end+1,1:2) = ls{1};
                    end
                end
            end
            fclose(fi);
        end %loadSpectraSuite

        function [d, xtype, ytype] = loadSpreadsheet(filename, args)
            % read table from file
            t = readtable(filename,args{:});
            d = t{:,:};
            if iscell(d)
                % convert to numeric array
                d = cellfun(@str2double,d);
                for dcol = 1:size(d,2)
                    % remove rows containing NaN
                    inan = isnan(d(:,dcol));
                    if any(inan)
                        d(inan,:) = [];
                    end
                end
            end
            xtype = t.Properties.VariableNames{1};
            if size(d,2)>2
                ytype = t.Properties.VariableNames(2:end);
            else
                ytype = t.Properties.VariableNames{2};
            end
        end %loadSpreadsheet
    end %private
end
