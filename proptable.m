mainprops = {'ID'; 'DateTime'; 'ExpID'};

% get property names of spectra
props = fields(SP(1));

% find properties that are not numeric arrays
pscalar = false(numel(props),1);
for k = 1:numel(props)
    prop = SP(1).(props{k});
    if isnumeric(prop) && ~isscalar(prop)
        pscalar(k) = false;
    else
        pscalar(k) = true;
    end
end
props = props(pscalar & ~ismember(props,mainprops));
props = [mainprops; props];

% Create a blank table
nprop = numel(props); % number of columns (properties)
nspec = numel(SP); % number of rows (spectra)
propcell = cell(nspec,nprop);

for k = 1:numel(SP)
    for p = 1:nprop
        propcell{k, p} = SP(k).(props{p});
    end
end

tbl = cell2table(propcell, 'VariableNames', props);