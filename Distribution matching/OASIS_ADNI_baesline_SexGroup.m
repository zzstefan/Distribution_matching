% In this script, we combine ADNI2 and OASIS-3 dataset with only baseline
% subjects, healthy subjects in OASIS3, healthy,EMCI,SMC subjects in ADNI2;
clc;
clear all;

dataTable = readtable('base_OASISADNI.csv', VariableNamingRule='preserve');

% Split data by sex
male_indices = strcmp(dataTable.Sex, 'M');
female_indices = strcmp(dataTable.Sex, 'F');

% Create separate tables for each sex
dataTable_male = dataTable(male_indices, :);
dataTable_female = dataTable(female_indices, :);

fprintf('Processing male group...\n');
process_sex_group(dataTable_male, 'male');

% Process female group
fprintf('Processing female group...\n');
process_sex_group(dataTable_female, 'female');

% Print summary of analysis
fprintf('\nAnalysis Summary:\n');
fprintf('Total subjects: %d\n', height(dataTable));
fprintf('Male subjects: %d\n', height(dataTable_male));
fprintf('Female subjects: %d\n', height(dataTable_female));

% For each site, show sex distribution
sites = unique(dataTable.SITE);
fprintf('\nSex distribution by site:\n');
for i = 1:length(sites)
    site = sites{i};
    site_indices = strcmp(dataTable.SITE, site);
    males = sum(site_indices & male_indices);
    females = sum(site_indices & female_indices);
    fprintf('%s: %d males, %d females\n', site, males, females);
end
