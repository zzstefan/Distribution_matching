% Script to combine and harmonize ADNI2, OASIS-3, and PREVENT-AD datasets
clc;
clear all;

% Load ADNI and OASIS data
dataTable = readtable('../DM_OASISADNI/base_OASISADNI.csv', VariableNamingRule='preserve');
Sc_ADNIOASIS = read_sc(dataTable, 'C');



% Load PREVENT-AD data
num_preventad = 340; % Replace with your actual indices
Sc_PREVENTAD = zeros(num_preventad, size(Sc_ADNIOASIS, 2));
load('prevent_AD_C.mat', 'c');
for i = 1:num_preventad
    % Extract raw connectivity matrix for subject i
    subject_matrix = c(:, :, i, 1);
    % Vectorize the connectivity matrix
    Sc_PREVENTAD(i,:) = vectorize_sc(subject_matrix);
end

% Combine Sc matrices and create site labels
Sc = [Sc_ADNIOASIS; Sc_PREVENTAD];
site_labels = [dataTable.SITE; repmat({'PREVENT-AD'}, num_preventad, 1)];

% Process combined data
Nonezero_col = find(any(Sc, 1));
sc = log(Sc(:,Nonezero_col)+1);
num_site = 3; % OASIS-3, ADNI2, and PREVENT-AD

% Initialize parameters
weight = zeros(size(sc,2), num_site);
params = zeros(size(sc,2), num_site, 2);

% Fit distributions for all features
for i=1:size(sc,2)
    [weight(i,:), params(i,:,:)] = fit_dist_new(sc, site_labels, i);
end

% Function to handle harmonization for different reference sites

reference_site = 'OASIS-3'; % 'PREVENT-AD' or 'OASIS-3'
harmonized_fea = harmonize_DM(sc, site_labels, params, weight, reference_site);

% Visualization
for i = 1:100
    fit_dist_histogram_threeDatasets(sc,harmonized_fea, site_labels, i, 0, reference_site);
end


har_fea = exp(harmonized_fea) - 1;
har_fea(isinf(har_fea)) = 0;

% Reconstruct full connectivity matrices
harmonized_data = zeros(size(sc));
harmonized_data(:,Nonezero_col) = har_fea;

% Save results
save_harmonized_results(harmonized_to_OASIS, 'OASIS', dataTable, preventad_indices);
save_harmonized_results(harmonized_to_PREVENTAD, 'PREVENTAD', dataTable, preventad_indices);
