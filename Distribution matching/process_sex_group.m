function process_sex_group(dataTable_group, sex_label)
    % Read and process structural connectivity matrices
    Sc = read_sc(dataTable_group, 'C');
    Nonezero_col = find(any(Sc, 1));
    sc = log(Sc(:,Nonezero_col)+1);

    num_site = length(unique(dataTable_group.SITE));
    weight = zeros(size(sc,2), num_site);
    params = zeros(size(sc,2), num_site, 2);

    % Fit gamma distribution
    for i = 1:size(sc,2)
        [weight(i,:), params(i,:,:)] = fit_dist(sc, dataTable_group, i);
    end

    % Perform harmonization
    harmonized_fea = harmonize_DM_ADNIOASIS(sc, dataTable_group, params, weight);

    % Create directory for harmonized matrices if it doesn't exist
    harmonized_dir = './harmonized_matrices_sex';
    if ~exist(harmonized_dir, 'dir')
        mkdir(harmonized_dir);
    end

    % Generate distribution plots with sex-specific folders
    for i = 1:100
        fit_dist_histogram_ADNIOASIS(sc, harmonized_fea, dataTable_group, i, 0, sex_label);
    end

    % Process harmonized features
    har_fea2 = harmonized_fea;
    har_fea = exp(har_fea2) - 1;
    har_fea(isinf(har_fea)) = 0;
    harmonized_data = zeros(size(Sc));
    harmonized_data(:,Nonezero_col) = har_fea;
    harmonized_matrix = vec_to_sym(harmonized_data, 85);

    % Save harmonized matrices with standardized naming
    for i = 1:height(dataTable_group)
        % Get subject ID
        subject_id = dataTable_group.Subject_ID{i};
        
        % Create standardized filenames
        c_filename = sprintf('%s_%s_DM_C.mat', subject_id, sex_label);
        ac_filename = sprintf('%s_%s_DM_aC.mat', subject_id, sex_label);
        
        % Save C matrix
        DM_C = squeeze(harmonized_matrix(i, :, :));
        save(fullfile(harmonized_dir, c_filename), 'DM_C');
        
        % Save augmented connectivity matrix
        DM_aC = augConnMatrix(DM_C);
        save(fullfile(harmonized_dir, ac_filename), 'DM_aC');
    end
end