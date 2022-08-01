function vars = extractvars(cindex)
% Extract variable structure from an index table
%
% vars = extractvars(cindex)
%
% Creates a struct vars with unique values found in the table cindex.
% 
% See also catindex

vars = table2struct(cindex,'ToScalar',true);
varnames = fields(vars);
for k = 1:numel(varnames)
    if iscategorical(vars.(varnames{k}))
        vars.(varnames{k}) = categories(vars.(varnames{k}));        
    else
        uvar = unique(vars.(varnames{k}));
        if(isnumeric(uvar))
            uvar(isnan(uvar)) = [];
        end
        vars.(varnames{k}) = uvar;
    end
end