function cindex = catindex(defs,vars,newvars)
% CATINDEX Create a categorical variable table based on list of 
% keywords and variable definitions or a "variable=value" definition list.
%
% Synthax
%     cindex = catindex(defs,vardef)
%     cindex = catindex(deflist)
%
% Description
% catindex(defs,vardef) searches for keywords in defs and sets the
% corresponding variable values as defined by vardef
%
% Categorical variables are dependent variables that define
% groups, such as 'Gender'. Variables have a
% defined set of values, such as 'Male', 'Female'.
% The struct vars contains fields for each variable.
% Each field must be a cell array of strings defining the
% possible values for the given variable.
% Alternatively, variables can be defined as a cell array of
% strings. Each string must contain a semicolon-separated values
%
% Input/Output parameters
%
%     defs - cell array of strings. Each string contains keywords that will
%            be assigned to variables
%     vardef - struct with fields that define variables. Each field
%            contains a number/cell array of possible values
%     cindex - a table with rows for each element in defs and columns
%              (variables) for each variable in vardef.
%
% catindex(deflist) does not need a separate variable definition but reads
%     variable names and corresponding values directly from deflist
%
%     deflist - cell array of strings with variable name/value definitions
%     
%     catindex parses the array of variable names and values and creates a
%     table with rows for each elements in the array as in the example.
%
% Examples
%
% vars.SampleGroup = { 'groupA', 'groupB', 'groupC' };
% vars.Treatment = {'control', 'treated'};
%
% cindex = catindex(filelist,vars);
% disp(cindex)
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
%
%     deflist = {'Name=Anna, Age=32, Gender=F';
%                'Name=John, Age=19, Gender=M'};
%     cindex = catindex(defs)
%     
% cindex = 
% 
%     Name    Age    Gender
%     ____    ___    ______
% 
%     Anna    32     F     
%     John    19     M      

if nargin == 1
    del1 = '=';
    del2 = ',';
    
    n = numel(defs);
    
    varstruct = struct.empty(n,0);
    
    rexp = strcat('([A-z]\w+)', del1, '([^', del1, del2, ']+)');
    rres = regexp(defs,rexp,'tokens');
    
    for k = 1:n
        if ~isempty(rres{k})
            rk = rres{k};
            for l = 1:numel(rk)
                name = rk{l}{1};
                value = rk{l}{2};
                numvalue = str2double(value);
                if strcmp(num2str(numvalue),value)
                    value = numvalue;
                end
                varstruct(k).(name) = value;
            end
        end
    end
    
    cindex = struct2table(varstruct);
    for v = 1:size(cindex,2)
        if iscellstr(cindex{:,v})
            cindex.(v) = categorical(cindex.(v));
        else
            str = arrayfun(@iscellstr,cindex{:,v});
            if any(str)
                cindex{~str,v} = {''};
                cindex.(v) = categorical(cindex.(v));
            end
        end
    end
elseif nargin >= 2
           cindex = struct;
           n = numel(defs);
           
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
               bindex.(varname) = repmat({defaultvar},n,1);    
               
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
                       if strfind(upper(defs{k}), upper(var))
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
           if numel(unique(defs)) == n
               cindex = struct2table(cindex,'RowNames',defs);
           else
               cindex = struct2table(cindex);
           end        
else
    error('Incorrect number of input arguments.')
end