function newdata = convert_specdata(olddata)
% convert_specdata(spdata) - update spdata from an older version 
%  of Spectr-O-Matic to the current
% 
% convert_specdata(filename) - update objects in .mat file from an older 
%  version of Spectr-O-Matic to the current

if ischar(olddata)
    counter = 0;
    if length(olddata>4) && ~strcmpi(olddata(end-3:end),'.mat')
        olddata = [olddata, '.mat'];
    end
    data = load(olddata);
    vars = fields(data);
    for k = 1:numel(vars)
        if isa(data.(vars{k}), 'SpecData')
            counter = counter + 1;
            data.(vars{k}) = convert_specdataobject(data.(vars{k}));
        end
    end
    if counter > 0
        bakname = [olddata(1:end-4), '.bak'];
        movefile(olddata, bakname);
        save(olddata, '-struct', 'data');
        fprintf('%d objects in %s updated. Original file saved as %s.\n', counter, olddata, bakname)
    else
        fprintf('No SpecData objects found in file %s.\n', olddata);
    end
else
    newdata = convert_specdataobject(olddata);
end

end

function newdata = convert_specdataobject(olddata)
props = {'ID','X','Y','XType','XUnit','YType','YUnit','ExpID','Comment','History','DateTime'};
n = numel(olddata);
newdata = specdata(n,0);
for k = 1:n
    for p = 1:numel(props)
        newdata(k).(props{p}) = olddata(k).(props{p});
    end
end
newdata = reshape(newdata, size(olddata));
end
    