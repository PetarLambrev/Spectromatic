classdef specdata2D < specparent
    % specdata2D - two-way (time-resolved) data
    % Spectr-O-Matic version 2.4
    % 
    % Petar Lambrev, 2012-2023
    %
    % Container for time-resolved spectroscopy data
    % A specdata2D dataset contains a 2D array of data (Z) dependent on two variables 
    % (X and T) plus additional identifying information, such as ID, units, etc. 
    % The class is intended for use with time-resolved spectroscopy data.
    %
    % X is an m*1 array containing the spectral axis (wavelength) values,
    % T is a 1*n array containing the time axis values, 
    % Y is an m*n array storing the measured quantity (abs, fluorescence, etc.). 
    %
    % specdata2D is derived on class SpecParent and shares common methods with
    % specdata (for steady-state spectroscopy data)
    %
    % See also: specdata

    properties
        T double = [];        % second axis (time)
        TType string = "T";   % time axis name
        TUnit string = "";    % time axis unit
    end
    
    properties (Dependent = true, SetAccess = private)
        dim
    end
    
    methods
       % Constructor
       function SP = specdata2D(newX,newT,newY,ID,varargin)
           % create specdata object
           %
           % Usage:
           % SP = specdata2D
           % Creates a blank spectrum SP
           %
           % SP = specdata2D(SParray)
           %
           % Convert a matrix of specdata
           % to a column vector of specdata2D
           %
           % SP = specdata2D(X, T, Y)
           % Creates a spectrum SP with supplied X, T and Y values
           %
           % SP = specdata2D(X, T, Y, ID)
           % Creates a spectrum SP with X, T, Y values and ID
           %
           % SP = specdata2D(X, T, Y, ID, ...)
           % Creates a spectrum SP with X, T, Y values, ID, and other
           % parameters, e.g. 'YType','CD','YUnit','mdeg',...
           %
           % The optional parameters can be also passed as one struct.
           
           proplist = {'ID','DateTime','ExpID','XType','XUnit','TType','TUnit','YType','YUnit','Comment','Metadata'};           
           
           if nargin==1 && isa(newX,'specdata')
              % Convert specdata to specdata2D
              nrows = size(newX,1);              
              ncols = size(newX,2);              
              
              SP = specdata2D.empty(nrows,0);
              for r = 1:nrows
                  % Copy properties from the source
                  newprops = properties(newX(r,1));
                  for p = 1:numel(proplist)
                      if ismember(proplist{p},newprops)
                          SP(r,1).(proplist{p}) = newX(r,1).(proplist{p});
                      end
                  end
                  SP(r,1).X = newX(r,1).X;
                  
                  % Determine number valid spectra in this row
                  dim = newX(r,1).dim;
                  validspec = false(1,ncols);
                  for c = 1:ncols
                      if newX(r,c).dim == dim
                          validspec(c) = true;
                      end
                  end
                  validspec = find(validspec);
                  
                  % Copy Y from all valid spectra into SP.Y
                  for y = 1:numel(validspec)
                      SP(r,1).Y(:,y) = newX(r,validspec(y)).Y;
                  end   
                  
                  % Create T array
                  SP(r,1).T = 1:numel(validspec);
              end
              % End constructor
           else   
               % Create specdata2D from arrays
               if (nargin < 3)
                   newX = [];
                   newT = [];
                   newY = [];
               end
           % validate X, T, Y
               if (~isreal(newX))
                   error('Cannot create specdata2D object. Invalid X data.')
               end
               if (~isreal(newT))
                   error('Cannot create specdata2D object. Invalid T data.')
               end
               if (~isreal(newY))
                   error('Cannot create specdata2D object. Invalid Y data.')
               end               
               if (length(newX)~=size(newY,1))
                   error('Cannot create specdata2D object. X and Y data dimensions do not match.');
               end
               if (length(newT)~=size(newY,2))
                   error('Cannot create specdata2D object. T and Y data dimensions do not match.');
               end
               
               SP.X(:,1) = newX(:);
               SP.T(1,:) = newT(:);
               [SP.X, Xind] = sort(SP.X);
               SP.Y = newY(Xind, :);
               
               if (nargin >= 4)
                   SP.ID = ID;
               end %if
               if nargin > 4
                   if isstruct(varargin(1))
                       optlist = varargin(1);
                   else
                       optlist = struct(varargin{:});
                   end
               else
                   optlist = [];
               end
               if (exist('optlist','var') && isstruct(optlist))
                   for i = 1:length(proplist)
                       if(isfield(optlist,proplist{i}))
                           SP.(proplist{i}) = optlist.(proplist{i});
                       end
                   end %for
               end %if               
           end %if
           
        end %Constructor

        %% Queries
       function n = get.dim(SP)
           % number of data points in X and T
           %
           % Usage:
           %    n = MySpectrum.dim
           n = zeros(length(SP),2);
           for i = 1:length(SP)
               n(i,1) = length(SP(i).X);
               n(i,2) = length(SP(i).T);
           end %for
       end %size
       
       function res = tind(SP,t)
           % TIND index of T values
           %
           % Synthax
           %    ti = A.tind(t)
           %
           % Returns the index of the T value(s) specified by t. 
           % If no exact match is found, returns the index of the closest T
           % If A is an array, ti will be a matrix
           % with rows corresponding to spectra and columns for each t
           res = zeros(numel(SP),numel(t));
           for k = 1:numel(SP)
               for l = 1:numel(t)
                   dif = abs(SP(k).T-t(l));
                   [~, res(k,l)] = min(dif(:));
               end
           end %for
       end %tind
       
       function res = Yt(SP,t)
           % YT Extract Y(T)     
           %
           % Synthax
           %    y = A.Yt(t)
           %
           % Description
           % Extracts the scalar values Y(t)
           % If A is an array, y is also an array 
           % containing Y values for each spectrum in A.
           %
           % Input parameters
           %    A - specdata2D object or array
           %    t - scalar value of T or array of scalar values
           %
           % See also: Yx
           Xn = size(SP(1).Y,1);
           res = zeros(Xn, length(SP));
           for k = 1:length(SP)
               if (size(SP(k).Y,1) ~= Xn)
                   error('The X data in the collection have different dimensions.');
               end
               Tind = SP(k).tind(t);
               if (isempty(Tind)), error('T value not found.'); end
               res(:,k) = SP(k).Y(:,Tind);
           end %for
       end %Yt    

       function res = Yxt(SP,x,t)
           % YT Extract Y(X,T)    
           %
           % Synthax
           %    y = A.Yxt(x,t)
           %
           % Description
           % Extracts the scalar values Y(x,t)
           % If A is an array, y is also an array 
           % containing Y values for each spectrum in A.
           %
           % Input/Output parameters
           %    A - specdata2D object or array. All spectra in A must have
           %    the same dimensions (X, T).
           %    x - scalar value of X or array of scalar values
           %    t - scalar value of T or array of scalar values       
           %
           %    y - a matrix of y values where each row contains the
           %        the Y(x,t) values for one spectrum in A
           %
           %    At least one of x or t must be a scalar.
           %

           % See also: Yx
           if isscalar(x)
               ncol = numel(t);
           elseif isscalar(t)
               ncol = numel(x);
           else
               error('At least one of x, t must be scalar.')
           end
           res = zeros(numel(SP),ncol);
           for k = 1:numel(SP)
               xi = SP(k).xind(x);
               ti = SP(k).tind(t);
               res(k,:) = SP(k).Y(xi,ti);
           end %for
       end %Ytx    
       
       %% Other calculations
       function res = settlim(SP,Tlim)
           % SETTLIM Set minimum and maximum limits for the T values.
           % 
           % Synthax
           %    NewData = MyData.settlim(Tlim)
           %
           % Description
           % The parameter Tlim must be an array of two elements: 
           % Tlim = [Tmin Tmax]. 
           % Any datapoints lying outside the specified range are deleted.
           res = SP;
           for i = 1:length(SP)
               j = 1;
               while (j <= length(res(i).T))
                   if (res(i).T(j) < Tlim(1)) || (res(i).T(j) > Tlim(2))
                       res(i).T(j) = [];
                       res(i).Y(:,j) = [];
                   else
                       j = j + 1;
                   end %if
               end %while
           end %for
       end %setXlim
       
function res = sett(SP,newT)
           % SETT Convert T axis
           %
           % Synthax
           %    NewSpectrum = MySpectrum.sett(newT)
           %
           % Replaces the original T axis values with newT and
           % redistributes the data. The Y values are interpolated 
           % to correspond to the new T array.
           % The sett function is useful for operations on spectra with 
           % different T scaling. To perform such operation, the T axes must be 
           % made identical first, by using an expression of the form
           %    NewSpectrum1 = MySpectrum1.setT(MySpectrum2.T)
           %
           % sett uses linear interpolation to match the Y values 
           % and does NOT extrapolate. If the new T array contains values outside
           % the range of the original T array, the corresponding Y values 
           % will be NaN. To avoid this, make sure the new T array does not 
           % exceed the range of the original T array. If necessary, 
           % use the setXlim method to limit the data range.
           %
           % See also: settlim
           classname = mfilename('class');
           if(isa(newT,classname))
               newT = newT.T;
           end
           res = SP;
           for i = 1:numel(SP)
               res(i).T = newT(:);
               res(i).Y = zeros(size(SP(i).Y,1),numel(newT));
               for j = 1:size(SP(i).Y,1)
                   res(i).Y(j,:) = interp1(SP(i).T,SP(i).Y(j,:),newT(:));
               end
           end %for
       end %sett       
       
       function res = swapxt(SP)
           % SWAPXT Exchange X and T axes (transpose Y)
           % 
           % Synthax
           %   B = A.swapxt
           %   B = swapxt(A)
           %
           % Description
           %  Exchanges the X and T axes of specdata2D spectra in A.           
           %  The Y values matrix is transposed accordingly.
           res = SP;
           for k = 1:numel(SP)
               res(k).X = SP(k).T;
               res(k).XType = SP(k).TType;
               res(k).XUnit = SP(k).TUnit;
               res(k).T = SP(k).X;
               res(k).TType = SP(k).XType;
               res(k).TUnit = SP(k).XUnit;
               res(k).Y = SP(k).Y';
           end
       end %swapxt
       
       function res = bint(SP,step)
           % Box-car aggregation across T data points
           %
           % Usage:
           %    NewData = MyData.bint(step)
           % 
           % Bins (aggregates) datapoints using the specified step. 
           % For example if MySpectrum has datapoints at every 0.1 ns,
           % MySpectrum.bint(10) will return a spectrum with datapoints every 1 ns,
           % where each new datapoint is an average of 10 original datapoints.
           % The bin function is useful to reduce the file size and noise 
           % in the data but it also decreases spectral resolution.
           %
           % See also: bin

           res = SP;           
           for i = 1:length(SP)
               Xn = length(SP(i).X);
               Tn = length(SP(i).T);
               bins = fix(Tn/step);
               res(i).T = zeros(bins,1);
               res(i).Y = zeros(Xn,bins);
               for j = 1:bins
                   ia = (j-1)*step + 1;
                   ib = j*step;
                   res(i).T(j) = mean(SP(i).T(ia:ib));
                   for k = 1:Xn
                       res(i).Y(k,j) = mean(SP(i).Y(k,ia:ib));
                   end
               end %for
               res(i).Comment = sprintf('%s; bint %.f', SP(i).Comment,step);
           end %for
       end %bin

       function T = table(SP, varargin)
           % TABLE Convert specdata to table
           %
           % Synthax
           %    T = table(S)
           %
           % Description
           % Creates a table with rows for every spectrum and columns
           % containing properties (metadata), such as ID, XType, etc.
           % The actual data (X, Y) are also included in the table.

           for ix = 1:length(SP)
               A = SP(ix);
               B = struct;
               B.ID = string(A.ID);
               B.DateTime = datetime(A.DateTime);
               B.ExpID = string(A.ExpID);
               B.dim = A.dim;
               B.X = A.X;
               B.Y = A.Y;
               B.T = A.T;
               B.XType = string(A.XType);
               B.XUnit = string(A.XUnit);               
               B.YType = string(A.YType);
               B.YUnit = string(A.YUnit);
               B.TType = string(A.TType);
               B.TUnit = string(A.TUnit);
               B.Comment = A.Comment;
               B.History = string(A.History);

               S(ix) = B;
           end
           T = struct2table(S);
       end
       
      %% Save
      function saveTIMP(SP,filename)
          savetimp(SP,filename)
      end
      
      function savetimp(SP,filename)
        % Save spectra as ASCII file formatted for TIMP (wavelength-explicit).
        %
        % Usage:
        %    dat.save(filename)
        % 
        % If dat contains more than one dataset, they are saved as separate files. 
        % filename can be a single string or a cell array of strings: 
        % {?file1?, ?file2?, ?} for each dataset. 
        % If a single file name is given, numbers will be appended to the filename 
        % to identify each dataset.
        %
        % To save the data as a MATLAB file (.mat) instead, 
        % use MATLAB's built-in save command:
        %    save filename dat 

          for i=1:length(SP)
              s = size(SP(i).Y,1);
              if(nargin>=2)
                  if(iscell(filename))
                    fout = filename{i};
                  else
                      if (i > 1)
                        fout = strcat(filename,sprintf('%g.ascii',i));
                      else
                        fout = filename;
                      end
                  end
              else
                  fout = strcat(SP(i).ID,'.ascii');
              end
              fprintf('Saving %s... ',fout);
              fid = fopen(fout,'wt');
              fprintf(fid,'%s\n',SP(i).ID);
              fprintf(fid,'%s\n',SP(i).Comment);
              fprintf(fid,'%s\n','Wavelength explicit');
              fprintf(fid,'Intervalnr %d\n\t',s);
              fprintf(fid,'%0g\t',SP(i).X(1:s-1));
              fprintf(fid,'%0g\n',SP(i).X(s));
              for j = 1:length(SP(i).T)
                  fprintf(fid,'%0g\t',SP(i).T(j));
                  fprintf(fid,'%0g\t',SP(i).Y(1:s-1,j));
                  fprintf(fid,'%0g\n',SP(i).Y(s,j));
              end
              fclose(fid);
              fprintf('Done.\n');
          end
      end
      
      function saveh5(SP, FileName)
           % SAVEH5 Save spectra as an HDF5 file.
           %
           % Synthax
           %    A.saveh5(filename)
           Profile = 'Simple';
           LDel = ' ['; RDel = ']';
           
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
           
           fid = H5F.create(FileName); 
           dtype = 'H5T_NATIVE_DOUBLE'; dcpl = 'H5P_DEFAULT';
           
           for k = 1:numel(SP)
              gid = H5G.create(fid,SP(k).ID,3); 
 
              X = SP(k).X;
              sid = H5S.create_simple(1, length(X), length(X));
              xid = H5D.create(gid,SP(k).XType,dtype,sid,dcpl);
              H5D.write(xid,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',dcpl,X);
              H5DS.set_label(xid,0,SP(k).XUnit);
              H5S.close(sid);                                          
              
              T = SP(k).T;
              sid = H5S.create_simple(1, length(T), length(T));
              tid = H5D.create(gid,SP(k).TType,dtype,sid,dcpl);
              H5D.write(tid,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',dcpl,T);
              H5DS.set_label(tid,0,SP(k).TUnit);
              H5S.close(sid);                                          

              Y = SP(k).Y;
              dims = [length(X), length(T)];
              sid = H5S.create_simple(2, dims, dims);
              yid = H5D.create(gid,SP(k).YType,dtype,sid,dcpl);
              H5D.write(yid,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',dcpl,Y');
              H5DS.set_label(yid,0,SP(k).XUnit);
              H5DS.set_label(yid,1,SP(k).TUnit);
              H5S.close(sid);                            
              
              H5DS.set_scale(xid,SP(k).XType);
              H5DS.attach_scale(yid,xid,0);
              H5DS.set_scale(tid,SP(k).TType);
              H5DS.attach_scale(yid,tid,1);              
              
              H5D.close(xid);
              H5D.close(tid);
              H5D.close(yid);
              H5G.close(gid);              
           end
%               h5create(FileName, [DS '/X'], [length(SP(k).X) 1])
%               h5write(FileName, [DS '/X'], SP(k).X(:))
%               h5create(FileName, [DS '/T'], [1 length(SP(k).T)])
%               h5write(FileName, [DS '/T'], SP(k).T(:)')
%               h5create(FileName, [DS '/Y'], fliplr(size(SP(k).Y)))
%               h5write(FileName, [DS '/Y'], SP(k).Y')
%               
%               fid = H5F.open(FileName);
%               did = H5D.open(fid, [DS '/X']);
%               H5DS.set_scale(did, 'Wavelength');
%               H5DS.set_label(did, [SP(k).XType LDel SP(k).XUnit RDel]);
%               did = H5D.open(fid, [DS '/T']);
%               H5DS.set_scale(did, 'T');
%               H5DS.set_label(did, [SP(k).TType LDel SP(k).TUnit RDel]);
%               did = H5D.open(fid, [DS '/Y']);
%               H5DS.set_label(did, [SP(k).YType LDel SP(k).YUnit RDel]);
%               
%               for p = 1:length(Props)                  
%                   h5writeatt(FileName, DS, Props{p}, SP(k).(Props{p}))
%               end
%            end
%            
%            h5writeatt(FileName,'/','Creator','spectromatic v2.3');
%            h5writeatt(FileName,'/','Profile',Profile)
%            
            H5F.close(fid)
            fprintf('File %s saved.\r',FileName);           
      end
      
      %% Plotting
      function plot2D(SP,viewmode)
          % plot data as a 2D colour map          
          if(length(SP)>10)
              warning('Only the first 10 datasets will be plotted.');
              SP = SP(1:10);
          end      

          if numel(SP)==1
              plotOneSP(SP)
          else
              for k = 1:numel(SP)
                  figure;
                  plotOneSP(SP(k))
                  title(SP(k).ID)
              end
          end
          
          function plotOneSP(SP1)
              % Cmax =  max(abs(SP.Y(:)))/1.5;
              % cl = length(colormap)/2;
              Z = SP1.Y';
              % C = fix((Z ./ Cmax + 1) .* cl);
              surf(SP1.X,SP1.T,Z); % ,C,'CDataMapping','direct');
              if (SP1.X(end)<SP1.X(1))
                  set(gca,'Xdir','rev');
                  xlim([SP1.X(end) SP1.X(1)]);
              else
                  xlim([SP1.X(1) SP1.X(end)]);
              end
              ylim([SP1.T(1) SP1.T(end)]);
              % set(gca,'YDir','rev');
              shading interp
              if (exist('viewmode','var'))
                  if viewmode == 3 || viewmode == 2
                      view(viewmode)
                  end
              else
                  view(2)
              end
              
              xlabel(sprintf('%s (%s)',SP1.XType, SP1.XUnit));
              ylabel(sprintf('%s (%s)',SP1.TType, SP1.TUnit));
              zlabel(sprintf('%s (%s)',SP1.YType, SP1.TUnit));
              colorbar
          end
      end %2D
      
      function plotxy(SP,SelT,varargin)
          % plot spectra (X vs Y) at selected T points
          %
          % Usage:
          %    dat.plotxy(selT,Options)
          %
          % Time-resolved spectra (Y = f(X)) are plotted at T = selT.
          % If exact match is not found, the next larger time point is plotted.
          %
          % Optional parameters:
          %    'SmoothLine' or 'SmoothLine','yes' - spectra are plotted as smooth lines
          %    'SmoothLine', 'no' - spectra are plotted as straight lines
          %
          % See also: plotTransients
          if (length(SP)>10)
              SP = SP(1:10);
              warning('Only the first 10 datasets will be plotted.');
          end
          
          if ~exist('SelT','var')
              SelT = [];
          end
          
          % Optional parameters 
           ShowExpID = false;
           NoTeX = false;
           SmoothL = 0;
           if nargin > 2
               for i = 1:nargin-2
                   switch(lower(varargin{i}))
                       case 'smoothline'
                           % if 'smooth' is the last argument then 
                           % assume smooth = yes
                           if (nargin>(i+2))
                               smootharg = lower(varargin{i+1});
                           else
                               smootharg = 'yes';
                           end
                           switch smootharg                               
                               case 'no'
                                   SmoothL = -1;
                               case 'yes'
                                   SmoothL = 1;
                           end
                   end %switch
               end %for
           end %if           
          
           if numel(SP)==1
               plotOneSP(SP)
           else
               for k = 1:numel(SP)
%                  figure;
                   plotOneSP(SP(k))
                   title(SP(k).ID)
               end
           end

          function plotOneSP(SP1)
              if ~isempty(SelT)
                  fy = specdata2D.ind(SP1.T, SelT);
              else
                  fy = 1:numel(SP1.T);
              end
              if(isempty(fy)), error 'The selected T values are not found in the dataset.'; end
              if SmoothL == 1 || (SmoothL == 0 && numel(SP1.X) < 24)
                  nx = (numel(SP1.X)-1) * 10;
                  dx = (max(SP1.X) - min(SP1.X)) / nx;
                  try
                      plX = min(SP1.X):dx:max(SP1.X);
                      plY = interp1(SP1.X,SP1.Y(:,fy),plX,'spline');
                      if numel(SP1.X) < 24
                          plot(plX,plY,'-o','MarkerIndices',1:10:numel(plX));
                      else
                          plot(plX,plY);
                      end
                  catch
                      warning('Cannot interpolate data for smoothing.')
                      plot(plX,plY);
                  end
              else
                  plot(SP1.X,SP1.Y(:,fy));
              end
              % create legends
              legends = cell(numel(fy,0));
              for li = 1:numel(fy)
                  legends{li} = sprintf('%0.2g %s', SP1.T(fy(li)), SP1.TUnit);
              end
              
              legend(legends,'location','best');
              if strcmp(SP1.XUnit,'')
                  Lbl = sprintf('%s',SP1.XType);
              else
                  Lbl = sprintf('%s (%s)',SP1.XType,SP1.XUnit);
              end
              xlabel(Lbl);
              if strcmp(SP1.YUnit,'')
                  Lbl = sprintf('%s',SP1.YType);
              else
                  Lbl = sprintf('%s (%s)',SP1.YType,SP1.YUnit);
              end
              ylabel(Lbl);
              Xmin = min(SP1.X); Xmax = max(SP1.X);
              xlim([Xmin Xmax])              
          end
      end %XY
      
      function plotty(SP,SelX)
          % plot time traces (T vs Y) at selected X values
          %
          % Usage:
          %    dat.plotty(X)
          %
          % If X is a vector time traces (Y = f(T)) are plotted for each X.
          % If exact X match is not found, the next larger X value in dat
          % is selected. If the function is called without parameters, 
          % all time traces in dat are plotted.
          %
          % See also: plotSpectra

          if (length(SP)>10)
              SP = SP(1:10);
              warning('Only the first 10 datasets will be plotted.');
          end

           if numel(SP)==1
               plotOneSP(SP)
           else
               for k = 1:numel(SP)
                   figure;
                   plotOneSP(SP(k))
                   title(SP(k).ID)
               end
           end          
          
          function plotOneSP(SP1)
              if(~exist('SelX','var'))
                  SelX = SP1.X;
              end
              for i1 = 1:length(SelX)
                  plot(SP1.T,SP1.Y(find(SP1.X >= SelX(i1),1),:));
                  hold all
                  l{i1} = sprintf('%g %s', SelX(i1),SP1.XUnit);
              end
              legend(l,'location','NorthEast','FontSize',11);
              xl = SP1.TType; if ~isempty(SP1.TUnit), xl = [xl, ' (',SP1.TUnit,')']; end
              xlabel(xl);
              yl = SP1.YType; if ~isempty(SP1.YUnit), yl = [yl, ' (',SP1.YUnit,')']; end
              ylabel(yl);
          end

      end %Transients

    end

    methods (Static)
        function spect = load(filepattern, options, firstrecord, lastrecord)
            % load data from ASCII text file(s)
            %
            % Usage:
            %    dat = specdata2D.load(filepattern, options)
            %
            % The files must contain only numerical data arranged in columns. 
            % The first column of the data must contain time points 
            % The following columns contain the measured data, 
            % where each column represents a wavelength channel.
            %
            % If time points occur several times (e.q. in repeated scan
            % measurements), the data are averaged.
            
            % filepattern ? can specify one or more filenames (e.g. ?*.DAT?). 
            %    If several files match the criteria, then dat is a cell
            %    array of datasets ? dat(1), dat(2), ? - for each file.
            % options ? a structure specifying optional parameters (as fields): 
            %    XType, XUnit, ExpID, Comment. 
            %    The X axis also can be specified here, 
            %    otherwise it will contain consecutive number 
            %    (X values are not read from the file). 
            %    IDs are not specified here but will carry the file names.
            
            if (iscell(filepattern))
                files = filepattern;
            else
                files = getdir(filepattern);
            end %if
            if (~iscell(files)) , error('No files selected.'); end
            if  (length(files)<1), error('No files selected.'); end
            fn = length(files);
            %data = cell(1,fn);
            for i1 = 1:fn
                fprintf('Loading file %s... ',files{i1});
                store = load(files{i1});
                % check number of time points
                timepoints = unique(store(:,1));
                tn = length(timepoints);
                cn = size(store,2); % number of columns
                rn1 = size(store,1)/tn; % number of records (scans)
                rn = fix(rn1);                
                if (rn1 > rn)
                    warning('Partial record will be ignored.');
                end
                if (rn < 1)
                    error('No valid records found in file.');
                end
                if (rn > 1)                 
                    % average if many records
                    jstart = 1;
                    jend = rn;
                    if (nargin>=3)
                        jstart = firstrecord;
                    end
                    if (nargin>=4)
                        jend = lastrecord;
                    end;                    
                    idata = zeros(tn,cn);
                    idata(:,1) = store(1:tn,1);
                    for j = jstart:jend
                        idata(:,2:cn) = idata(:,2:cn) + store(j*tn-tn+1:j*tn,2:cn);
                    end;
                    idata(:,2:cn) = idata(:,2:cn)/rn;
                    data = idata;
                else
                    % only one record
                    data = store;                    
                end;
                fprintf('Done, %g record(s) read, %g time points.\n',rn, tn);
                cols = size(data,2);
                T = data(:,1)';
                Y = data(:,2:cols)';
                X = (1:cols-1)';
                spect(i1) = specdata2D(X,T,Y,files{i1});
                if (nargin < 2)
                    options = struct;
                end
                props = {'X','XType','XUnit','TType','TUnit','YType','YUnit','ExpID','Comment','Date'};
                if (exist('options','var'))
                    for p = 1:length(props)
                        if (isfield(options,props{p}))
                            spect(i1).(props{p}) = options.(props{p});
                        end
                    end %for p
                end %if exist
            end
        end %loadData
        
    end %Static
    
    %% Private methods
    methods (Access = protected, Hidden = true)
        function match = comp(SP, SP2)
            % compare two spectra if their axes match
            % return 0 if there is no difference, 1 if there is a
            % difference in X, 2 if there is a difference in Y
            match = 0;
            if ~isequal(SP.X,SP2.X) 
                match = 1;
            end
            if ~isequal(SP.T,SP2.T)
                match = 2;
            end            
        end        
    end % private methods

end % classdef
    