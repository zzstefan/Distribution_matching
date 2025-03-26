% In this script, we combine ADNI2 and OASIS-3 dataset with only baseline
% subjects, healthy subjects in OASIS3, healthy,EMCI,SMC subjects in ADNI2;


dataTable = readtable('base_OASISADNI.csv', VariableNamingRule='preserve');

% delete those AD and LMCI from ADNI2 dataset;


%dataTable.Properties.VariableNames{'Site'}='SITE';
Sc = read_sc(dataTable,'C');
Nonezero_col = find(any(Sc, 1));
sc = log(Sc(:,Nonezero_col)+1);

num_site = length(unique(dataTable.SITE));
weight = zeros(size(sc,2), num_site);
params = zeros(size(sc,2), num_site, 2);


%fit the gamma distribution and estimate shape and scale parameters;
for i=1:size(sc,2)
    [weight(i,:), params(i,:,:)] = fit_dist(sc, dataTable, i);
end


%% distribution matching to harmonized features from ADNI2 to OASIS-3;
harmonized_fea = harmonize_DM_ADNIOASIS(sc,dataTable,params,weight);


% Save the results for correlation analysis;
%save('./Result/correlation_with_variables/DM_C.mat', 'harmonized_fea');

% save to fitted_distribution_all_matlab
% two columns, first is the unharmonized version and second is the harmonized
% version;
for i = 1:100
    fit_dist_histogram_ADNIOASIS(sc,harmonized_fea, dataTable, i, 0);
end


%% examine the result with the inconsistency in raw sc and dm sc;
% to_draw = load('/autofs/space/genesis_001/users/Projects/PPMI/code/Result/to_matlab/to_draw_DM.mat');
% for i = 1:10
%     fit_dist_histogram_ADNIOASIS(to_draw.raw_sc,to_draw.DM_sc, dataTable, i, 0);
% end


%% save the harmonized features to specific data folders
har_fea2 = harmonized_fea;
har_fea = exp(har_fea2) - 1; % log(unhar_fea+1)
har_fea(isinf(har_fea)) = 0;
harmonized_data = zeros(size(Sc));
harmonized_data(:,Nonezero_col) = har_fea;
harmonized_matrix = vec_to_sym(harmonized_data, 85);



%%%
% for i=1:size(harmonized_matrix,1)
%     aC(i,:,:) = augConnMatrix(squeeze(harmonized_matrix(i,:,:)));
% end
%%%

% Save each C matrix to a .mat file
dataTable.DM_C = strrep(dataTable.conn_mat, 'connMat_Hough.mat', 'baseline_ADNIOASIS/DM_C.mat');

for i = 1:height(dataTable)
    [folderPath, ~, ~] = fileparts(dataTable.DM_C{i});
    if ~exist(folderPath, 'dir')
    % If the folder does not exist, create it
        mkdir(folderPath);
    end
    DM_C = squeeze(harmonized_matrix(i, :, :));
    save(dataTable.DM_C{i}, 'DM_C');
end




%calcualte and save aC;
dataTable.DM_aC = strrep(dataTable.conn_mat, 'connMat_Hough.mat', 'baseline_ADNIOASIS/DM_aC.mat');
for i= 1:height(dataTable)
    DM_aC= augConnMatrix(squeeze(harmonized_matrix(i,:,:)));
    save(dataTable.DM_aC{i}, 'DM_aC');
end
