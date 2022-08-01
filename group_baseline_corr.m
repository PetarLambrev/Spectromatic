function corrdata = group_baseline_corr(data,id_group,id_baseline)

% group_baseline_corr - Subtract baselines from groups of spectra
% 
% corr = group_baseline_corr(data,id_group,id_baseline)
%
% data - SpecData object containing groups of spectra and baselines
% id_group - cell array of strings defining groups
% id_baseline - string or cell array of strings defining baseline groups
%
% For each group of spectra defined by the keyword id_group,
% there must be a corresponding baseline. If id_baseline is a cell array,
% then it must contain the same number of elements as id_group.
%
% Alternatively id_baseline can be a single string, then baselines are
% defined by id_baseline + id_group.

if size(data,1) == 1
    data = data';
end

for k = id_group
    g = data(:,1).fi(k{1});
    bl = g & data(:,1).fi(id_baseline);
    data(g,:) = data(g,:) - data(bl,:);
    data(bl,:) = [];
end
corrdata = data;