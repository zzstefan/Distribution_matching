% Function to process each sex group
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

    % Create output directory for sex-specific results
    result_dir = sprintf('./Result_%s', sex_label);
    if ~exist(result_dir, 'dir')
        mkdir(result_dir);
    end

    % Generate distribution plots
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

    % Save harmonized matrices
    for i = 1:height(dataTable_group)
        % Create sex-specific output path
        output_path = strrep(dataTable_group.conn_mat{i}, ...
            'connMat_Hough.mat', ...
            sprintf('baseline_ADNIOASIS_%s/DM_C.mat', sex_label));
        
        [folderPath, ~, ~] = fileparts(output_path);
        if ~exist(folderPath, 'dir')
            mkdir(folderPath);
        end
        
        % Save C matrix
        DM_C = squeeze(harmonized_matrix(i, :, :));
        save(output_path, 'DM_C');
        
        % Save augmented connectivity matrix
        aC_path = strrep(output_path, 'DM_C.mat', 'DM_aC.mat');
        DM_aC = augConnMatrix(DM_C);
        save(aC_path, 'DM_aC');
    end
    
    % Save summary statistics
    summary_stats = struct();
    summary_stats.num_subjects = height(dataTable_group);
    summary_stats.sites = unique(dataTable_group.SITE);
    summary_stats.subjects_per_site = histcounts(categorical(dataTable_group.SITE));
    summary_stats.weight = weight;
    summary_stats.params = params;
    
    %save(fullfile(result_dir, 'summary_stats.mat'), 'summary_stats');
end