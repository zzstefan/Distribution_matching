% Script to perform comprehensive experiments with ADNI, OASIS, and PREVENT-AD
% Four experiments:
% 1. No sex grouping, OASIS reference
% 2. No sex grouping, PREVENT-AD reference
% 3. Sex-specific grouping, OASIS reference
% 4. Sex-specific grouping, PREVENT-AD reference



% Create directory for results
main_output_dir = './harmonization_results';
if ~exist(main_output_dir, 'dir')
    mkdir(main_output_dir);
end

% Load ADNI and OASIS data
dataTable = readtable('../DM_OASISADNI/base_OASISADNI.csv', VariableNamingRule='preserve');
Sc_ADNIOASIS = read_sc(dataTable, 'C');

% Load PREVENT-AD data
num_preventad = 340;
Sc_PREVENTAD = zeros(num_preventad, size(Sc_ADNIOASIS, 2));
load('prevent_AD_C.mat', 'c');

% Process PREVENT-AD data
for i = 1:num_preventad
    subject_matrix = c(:, :, i, 1);
    Sc_PREVENTAD(i,:) = vectorize_sc(subject_matrix);
end

preventad_info = readtable('prevent_AD.csv', VariableNamingRule='preserve');

% Extract subject IDs and sex information
preventad_ids = preventad_info.Subject_ID;
preventad_sex = preventad_info.Sex;

% Define experiment parameters
reference_sites = {'OASIS-3', 'PREVENT-AD'};
sex_groups = [false, true]; % false = no sex grouping, true = sex grouping

% -------------------------------------------------------------------------
% EXPERIMENT PIPELINE FUNCTION
% -------------------------------------------------------------------------
function run_experiment(Sc_ADNIOASIS, Sc_PREVENTAD, dataTable, preventad_ids, preventad_sex, reference_site, sex_grouping, output_dir)
    fprintf('Running experiment: Reference=%s, Sex Grouping=%d\n', reference_site, sex_grouping);
    
    % Combine data from all three datasets
    Sc = [Sc_ADNIOASIS; Sc_PREVENTAD];
    site_labels = [dataTable.SITE; repmat({'PREVENT-AD'}, length(preventad_ids), 1)];
    
    % Process combined data
    Nonezero_col = find(any(Sc, 1));
    sc = log(Sc(:,Nonezero_col)+1);
    
    % Create experiment output directory
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    if ~sex_grouping
        % NO SEX GROUPING - PROCESS ALL DATA TOGETHER
        % -----------------------------------------
        
        % Fit distributions for all features
        num_site = 3; % OASIS-3, ADNI2, and PREVENT-AD
        weight = zeros(size(sc,2), num_site);
        params = zeros(size(sc,2), num_site, 2);
        
        for i = 1:size(sc,2)
            [weight(i,:), params(i,:,:)] = fit_dist_new(sc, site_labels, i);
        end
        
        % Perform harmonization
        harmonized_fea = harmonize_DM(sc, site_labels, params, weight, reference_site);
        
        % Visualization
        for i = 1:100
            fit_dist_histogram_threeDatasets(sc, harmonized_fea, site_labels, i, 0, reference_site);
        end
        
        % Convert back from log space
        har_fea = exp(harmonized_fea) - 1;
        har_fea(isinf(har_fea)) = 0;
        
        % Reconstruct full connectivity matrices
        harmonized_data = zeros(size(Sc));
        harmonized_data(:,Nonezero_col) = har_fea;
        
        % Save ADNI/OASIS results
        for i = 1:height(dataTable)
            subject_id = dataTable.Subject_ID{i};
            
            % Create output filename
            c_filename = sprintf('%s_DM_C.mat', subject_id);
            ac_filename = sprintf('%s_DM_aC.mat', subject_id);
            
            % Convert to full matrix and save
            DM_C = vec_to_sym(harmonized_data(i,:), 85);
            if ndims(DM_C) == 3
                % If the result is a 3D array of size 1×85×85, extract the 2D matrix
                DM_C = squeeze(DM_C);
            end
            save(fullfile(output_dir, c_filename), 'DM_C');
            
            % Save augmented matrix
            DM_aC = augConnMatrix(DM_C);
            save(fullfile(output_dir, ac_filename), 'DM_aC');
        end
        
        % Save PREVENT-AD results
        for i = 1:length(preventad_ids)
            idx = height(dataTable) + i;
            subject_id = preventad_ids{i};
            
            % Create output filename
            c_filename = sprintf('%s_DM_C.mat', subject_id);
            ac_filename = sprintf('%s_DM_aC.mat', subject_id);
            
            % Convert to full matrix and save
            DM_C = vec_to_sym(harmonized_data(idx,:), 85);
            if ndims(DM_C) == 3
                % If the result is a 3D array of size 1×85×85, extract the 2D matrix
                DM_C = squeeze(DM_C);
            end
            save(fullfile(output_dir, c_filename), 'DM_C');
            % Save augmented matrix
            DM_aC = augConnMatrix(DM_C);
            save(fullfile(output_dir, ac_filename), 'DM_aC');
        end
        
    else
        % SEX-SPECIFIC GROUPING - PROCESS MALE AND FEMALE SEPARATELY
        % -----------------------------------------
        
        % Combine sex labels
        all_sex_labels = [dataTable.Sex; preventad_sex];
        
        % Process male and female data separately
        for sex_val = {'M', 'F'}
            sex_label = sex_val{1};
            fprintf('  Processing %s group\n', sex_label);
            
            % Get sex-specific indices
            sex_indices = strcmp(all_sex_labels, sex_label);
            
            % Extract sex-specific data
            sc_sex = sc(sex_indices, :);
            site_labels_sex = site_labels(sex_indices);
            
            % Create sex-specific output directory
            sex_output_dir = fullfile(output_dir, lower(sex_label));
            if ~exist(sex_output_dir, 'dir')
                mkdir(sex_output_dir);
            end
            
            % Fit distributions for all features
            num_site = 3; % OASIS-3, ADNI2, and PREVENT-AD
            weight = zeros(size(sc_sex,2), num_site);
            params = zeros(size(sc_sex,2), num_site, 2);
            
            for i = 1:size(sc_sex,2)
                [weight(i,:), params(i,:,:)] = fit_dist_new(sc_sex, site_labels_sex, i);
            end
            
            % Perform harmonization
            harmonized_fea_sex = harmonize_DM(sc_sex, site_labels_sex, params, weight, reference_site);
            
            % Visualization
            for i = 1:100
                folder_name = sprintf('./fitted_histogram_to%s_%s/', reference_site, lower(sex_label));
                if ~exist(folder_name, 'dir')
                    mkdir(folder_name);
                end
                fit_dist_histogram_threeDatasets(sc_sex, harmonized_fea_sex, site_labels_sex, i, 0, reference_site, sex_label);
            end
            
            % Convert back from log space
            har_fea = exp(harmonized_fea_sex) - 1;
            har_fea(isinf(har_fea)) = 0;
          

            % Reconstruct full connectivity matrices - THIS WAS MISSING
            harmonized_data_sex = zeros(size(sc_sex, 1), size(Sc, 2));
            harmonized_data_sex(:, Nonezero_col) = har_fea;


            % Get subject IDs for this sex group
            all_subject_ids = [dataTable.Subject_ID; preventad_ids];
            subject_ids_sex = all_subject_ids(sex_indices);
            
            % Save harmonized matrices
            for i = 1:length(subject_ids_sex)
                subject_id = subject_ids_sex{i};
                
                % Create output filename
                c_filename = sprintf('%s_%s_DM_C.mat', subject_id, lower(sex_label));
                ac_filename = sprintf('%s_%s_DM_aC.mat', subject_id, lower(sex_label));
                
                % Convert to full matrix and save
                DM_C = vec_to_sym(harmonized_data_sex(i,:), 85);
                if ndims(DM_C) == 3
                    % If the result is a 3D array of size 1×85×85, extract the 2D matrix
                    DM_C = squeeze(DM_C);
                end
                save(fullfile(sex_output_dir, c_filename), 'DM_C');
                % Save augmented matrix
                DM_aC = augConnMatrix(DM_C);
                save(fullfile(sex_output_dir, ac_filename), 'DM_aC');
            end
        end
    end
    
    fprintf('Experiment completed: Reference=%s, Sex Grouping=%d\n\n', reference_site, sex_grouping);
end

% -------------------------------------------------------------------------
% RUN ALL EXPERIMENTS
% -------------------------------------------------------------------------

% Experiment 1: No sex grouping, OASIS reference
output_dir1 = fullfile(main_output_dir, 'no_sex_OASIS_ref');
run_experiment(Sc_ADNIOASIS, Sc_PREVENTAD, dataTable, preventad_ids, preventad_sex, 'OASIS-3', false, output_dir1);

% Experiment 2: No sex grouping, PREVENT-AD reference
output_dir2 = fullfile(main_output_dir, 'no_sex_PREVENTAD_ref');
run_experiment(Sc_ADNIOASIS, Sc_PREVENTAD, dataTable, preventad_ids, preventad_sex, 'PREVENT-AD', false, output_dir2);

% Experiment 3: Sex-specific grouping, OASIS reference
output_dir3 = fullfile(main_output_dir, 'sex_specific_OASIS_ref');
run_experiment(Sc_ADNIOASIS, Sc_PREVENTAD, dataTable, preventad_ids, preventad_sex, 'OASIS-3', true, output_dir3);

% Experiment 4: Sex-specific grouping, PREVENT-AD reference
output_dir4 = fullfile(main_output_dir, 'sex_specific_PREVENTAD_ref');
run_experiment(Sc_ADNIOASIS, Sc_PREVENTAD, dataTable, preventad_ids, preventad_sex, 'PREVENT-AD', true, output_dir4);

% -------------------------------------------------------------------------
% Helper function to vectorize connectivity matrix
% -------------------------------------------------------------------------
function vec = vectorize_sc(matrix)
    triu_idx = triu(true(size(matrix)), 1);
    vec = matrix(triu_idx)';
end