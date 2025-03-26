import numpy as np
import pandas as pd
import glob
import re
import scipy.stats as stats
from utils import *
import scipy.stats
from datetime import datetime
import scipy.io as sio
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import wilcoxon
from covbat_revise3 import apply_covbat
import os

# Load the datasets
# ADNI and OASIS in one file, Prevent-AD in another
adni_oasis_data = pd.read_csv('./matlab_code/DM_OASISADNI/base_OASISADNI.csv')
prevent_ad_data = pd.read_csv('./matlab_code/Prevent_AD/prevent_AD.csv')
prevent_ad_data['SITE'] = 'PREVENT-AD'

# Read the Prevent-AD data once from the 4D matrix
preventad_conn_matrix = scipy.io.loadmat('./matlab_code/Prevent_AD/prevent_AD_C.mat')['c']

# Combine the datasets
data = pd.concat([adni_oasis_data, prevent_ad_data], ignore_index=True)


def read_sc(data_df, matrix_type, preventad_matrix=None):
    matrices = np.zeros((len(data_df), int(85 * 84 / 2)))
    
    # Process ADNI and OASIS subjects
    adni_oasis_idx = data_df[data_df['SITE'].isin(['ADNI2', 'OASIS-3'])].index
    for i in adni_oasis_idx:
        try:
            tmp = scipy.io.loadmat(data_df.loc[i, 'conn_mat'])[matrix_type]
            matrices[i, :] = tmp[np.triu_indices(85, k=1)]
        except Exception as e:
            print(f"Error reading {matrix_type} for subject {data_df.loc[i, 'Subject_ID']}: {str(e)}")
    
    # Process Prevent-AD subjects - read directly from the 4D matrix
    preventad_idx = data_df[data_df['SITE'] == 'PREVENT-AD'].index
    dim_idx = 0 if matrix_type == 'C' else 1  # Index for C or aC in 4th dimension
    
    for i, idx in enumerate(preventad_idx):
        try:
            # Get matrix for this Prevent-AD subject from the 4D array
            tmp = preventad_matrix[:, :, i, dim_idx]
            matrices[idx, :] = tmp[np.triu_indices(85, k=1)]
        except Exception as e:
            print(f"Error extracting {matrix_type} for Prevent-AD subject at index {i}: {str(e)}")
    
    return matrices



# Create paths for distribution matching results based on your folder structure
def create_dm_paths(data_df, base_dir, ref_site, include_gender=False):
    """Create paths for distribution matching results with your folder structure"""
    paths = []
    
    if include_gender:
        # For sex-specific harmonization
        folder_name = f"sex_specific_{ref_site}_ref"
        for i, row in data_df.iterrows():
            sub_id = row['Subject_ID']
            sex = 'f' if row['Sex'].lower() == 'f' else 'm'  # Convert to lowercase 'f' or 'm'
            # Format: sex_specific_OASIS_ref/f/subjectID_f_DM_C.mat
            path = os.path.join(base_dir, folder_name, sex, f"{sub_id}_{sex}_DM_C.mat")
            paths.append(path)
    else:
        # For no-sex harmonization
        folder_name = f"no_sex_{ref_site}_ref"
        for i, row in data_df.iterrows():
            sub_id = row['Subject_ID']
            path = os.path.join(base_dir, folder_name, f"{sub_id}_DM_C.mat")
            paths.append(path)
            
    return paths




# Dictionary to store all results
results = {}

# First, run the combat and covbat harmonization once (no reference site needed)
# Identify non-zero columns for log transformation
C = read_sc(data, 'C', preventad_conn_matrix)  # Original connectivity
aC = read_sc(data, 'aC',preventad_conn_matrix)  # Original adjusted connectivity


# Identify non-zero columns for log transformation
nonZero_col = np.where(C.any(axis=0))[0]
zero_col = np.where(~C.any(axis=0))[0]
sc = np.delete(C, zero_col, axis=1)
sc2 = np.log(sc+1)

batch = data['SITE'].values
age = data['Age'].values
sex = data['Sex'].values

# Apply Combat harmonization methods
# 1. Combat with age
model, combat_age_C = harmonize(sc2, data, 0)
com_C = recon_full(combat_age_C, data, nonZero_col)

# 2. Combat without age
model, combat_noAge_C = harmonize_noAge(sc2, data)
com_noAge_C = recon_full(combat_noAge_C, data, nonZero_col)

# 3. Combat without age or sex
model, combat_noAge_noSex_C = harmonize_noAge_noSex(sc2, data)
com_noAge_noSex_C = recon_full(combat_noAge_noSex_C, data, nonZero_col)

# Same process for adjusted connectivity (aC)
# Combat harmonization with age as input on aC
model, com_aC = harmonize(np.log(aC+1), data, 0)

# Combat harmonization without age as input on aC
model, com_noAge_aC = harmonize_noAge(np.log(aC+1), data)

# Combat harmonization without age or sex as input on aC
model, com_noAge_noSex_aC = harmonize_noAge_noSex(np.log(aC+1), data)

# Apply CovBat harmonization methods
# For C matrices
cov_C_tmp = apply_covbat(sc2, batch, age, sex)
cov_C = recon_full(cov_C_tmp, data, nonZero_col)

cov_noAge_C_tmp = apply_covbat(sc2, batch, sex=sex)
cov_noAge_C = recon_full(cov_noAge_C_tmp, data, nonZero_col)

cov_noAge_noSex_C_tmp = apply_covbat(sc2, batch)
cov_noAge_noSex_C = recon_full(cov_noAge_noSex_C_tmp, data, nonZero_col)

# For aC matrices
cov_aC = apply_covbat(np.log(aC+1), batch, age, sex)
cov_noAge_aC = apply_covbat(np.log(aC+1), batch, sex=sex)
cov_noAge_noSex_aC = apply_covbat(np.log(aC+1), batch)

# Create a base SC dictionary with Combat and CovBat results (no reference site specific)
base_SC = {
    'C': norm_fea(C),
    'Com_C': norm_fea(np.exp(com_C)-1),
    'Com_noAge_C': norm_fea(np.exp(com_noAge_C)-1),
    'Com_noAge_noSex_C': norm_fea(np.exp(com_noAge_noSex_C)-1),
    'Cov_C': norm_fea(np.exp(cov_C)-1),
    'Cov_noAge_C': norm_fea(np.exp(cov_noAge_C)-1),
    'Cov_noAge_noSex_C': norm_fea(np.exp(cov_noAge_noSex_C)-1),
    'aC': norm_fea(aC),
    'Com_aC': norm_fea(np.exp(com_aC)-1),
    'Com_noAge_aC': norm_fea(np.exp(com_noAge_aC)-1),
    'Com_noAge_noSex_aC': norm_fea(np.exp(com_noAge_noSex_aC)-1),
    'Cov_aC': norm_fea(np.exp(cov_aC)-1),
    'Cov_noAge_aC': norm_fea(np.exp(cov_noAge_aC)-1),
    'Cov_noAge_noSex_aC': norm_fea(np.exp(cov_noAge_noSex_aC)-1)
}


# Calculate the base correlations with age
age_index = data.loc[~data.Age.isna()].index
age_values = data.loc[age_index, 'Age'].values

base_r = np.zeros((aC.shape[1], len(base_SC)))
base_p = np.zeros((aC.shape[1], len(base_SC)))

for j, c in enumerate(list(base_SC.keys())):
    fea = base_SC[c][age_index, :]
    for i in range(aC.shape[1]):
        base_r[i, j], base_p[i, j] = scipy.stats.pearsonr(age_values, fea[:, i])

# Apply Bonferroni correction
base_p = base_p * 3570

# Calculate and print statistics for base harmonization methods
R = np.abs(base_r)
mean_r = np.nanmean(R, axis=0)
std_r = np.nanstd(R, axis=0)

print("\n=== Age Correlation Results for Base Harmonization Methods ===")
print("Format: Mean (Std) of absolute correlations")
for i, key in enumerate(base_SC.keys()):
    print(f"{key}: {mean_r[i]:.2f} ({std_r[i]:.2f})")

# Wilcoxon tests comparing original vs harmonized
print("\n=== Wilcoxon Signed-Rank Tests for Base Harmonization Methods ===")
print("Comparing absolute correlations (|r|) between feature sets")

# Compare each harmonized version to original
base_indices = [0, 7]  # Indices for 'C' and 'aC'
for base_idx in base_indices:
    base_name = list(base_SC.keys())[base_idx]
    print(f"\nBase: {base_name}")
    for j in range(len(base_SC)):
        if j != base_idx and (j < 7 if base_idx == 0 else j >= 7):
            comp_name = list(base_SC.keys())[j]
            res = wilcoxon(R[:, j], R[:, base_idx], alternative='greater', nan_policy='omit')
            print(f"{comp_name} vs {base_name}: statistic={res.statistic:.2f}, p-value={res.pvalue:.0e}")

# Count significant correlations
sig_counts = (base_p < 0.05).sum(axis=0)
print("\n=== Number of Significant Age Correlations (p < 0.05) for Base Harmonization Methods ===")
for i, key in enumerate(base_SC.keys()):
    print(f"{key}: {sig_counts[i]}")

# Minimum p-values
min_p = np.nanmin(base_p, axis=0)
print("\n=== Minimum Bonferroni-Corrected p-values for Base Harmonization Methods ===")
for i, key in enumerate(base_SC.keys()):
    print(f"{key}: {min_p[i]:.0e}")


# Helper function to read DM matrices with proper handling for different datasets
def read_dm_matrices(data_df, dm_path_column, matrix_name):
    """Read distribution matching matrices for all datasets consistently"""
    matrices = np.zeros((len(data_df), int(85 * 84 / 2)))
    
    for i in range(len(data_df)):
        try:
            path = data_df.loc[i, dm_path_column]
            tmp = scipy.io.loadmat(path)[matrix_name]
            matrices[i, :] = tmp[np.triu_indices(85, k=1)]
        except Exception as e:
            print(f"Error reading {matrix_name} from {path} for subject {data_df.loc[i, 'Subject_ID']}: {str(e)}")
    
    return matrices

# Now run the distribution matching analysis for each reference site
# Path to the harmonization results
harmonization_dir = './matlab_code/Prevent_AD/harmonization_results'
# Define reference sites for harmonization
reference_sites = ['OASIS','PREVENTAD']
gender_options = [False,True]  # With or without gender grouping



for ref_site in reference_sites:
    for use_gender in gender_options:
        # Set up experiment identifier
        exp_id = f"DM_{ref_site}_{'withGender' if use_gender else 'noGender'}"
        
        # Create paths for distribution matching results
        data[f'DM_C_{exp_id}'] = create_dm_paths(data, harmonization_dir, ref_site, use_gender)
        data[f'DM_aC_{exp_id}'] = [path.replace('_DM_C.mat', '_DM_aC.mat') for path in data[f'DM_C_{exp_id}']]
        
        # Read distribution matching results with the new helper function
        DM_C = read_dm_matrices(data, f'DM_C_{exp_id}', 'DM_C')
        DM_aC = read_dm_matrices(data, f'DM_aC_{exp_id}', 'DM_aC')
        
        # Create SC dictionary with DM results for this reference site
        SC_dm = {
            'C': base_SC['C'],
            f'DM_C_{exp_id}': norm_fea(DM_C),
            'aC': base_SC['aC'],
            f'DM_aC_{exp_id}': norm_fea(DM_aC)
        }
        
        # Calculate correlation with age for DM methods
        r_dm = np.zeros((aC.shape[1], len(SC_dm)))
        p_dm = np.zeros((aC.shape[1], len(SC_dm)))
        
        for j, c in enumerate(list(SC_dm.keys())):
            fea = SC_dm[c][age_index, :]
            for i in range(aC.shape[1]):
                r_dm[i, j], p_dm[i, j] = scipy.stats.pearsonr(age_values, fea[:, i])
        
        # Apply Bonferroni correction
        p_dm = p_dm * 3570
        
        # Replace NaNs with 1.0 for p-values
        if np.isnan(p_dm).any():
            print(f"Warning: Found {np.isnan(p_dm).sum()} NaN p-values in {exp_id}. Replacing with 1.0")
            p_dm = np.nan_to_num(p_dm, nan=1.0)
        
        # Replace NaNs with 0.0 for r-values
        if np.isnan(r_dm).any():
            print(f"Warning: Found {np.isnan(r_dm).sum()} NaN r-values in {exp_id}. Replacing with 0.0")
            r_dm = np.nan_to_num(r_dm, nan=0.0)
        
        # Store results for this experiment
        results[exp_id] = {
            'r': r_dm,
            'p': p_dm,
            'SC_keys': list(SC_dm.keys())
        }
        
        # Calculate and print statistics for DM methods
        R_dm = np.abs(r_dm)
        mean_r_dm = np.nanmean(R_dm, axis=0)
        std_r_dm = np.nanstd(R_dm, axis=0)
        
        print(f"\n=== Age Correlation Results for {exp_id} ===")
        print("Format: Mean (Std) of absolute correlations")
        for i, key in enumerate(SC_dm.keys()):
            print(f"{key}: {mean_r_dm[i]:.2f} ({std_r_dm[i]:.2f})")
        
        # Wilcoxon tests comparing original vs harmonized for DM
        print(f"\n=== Wilcoxon Signed-Rank Tests for {exp_id} ===")
        print("Comparing absolute correlations (|r|) between feature sets")
        
        # Compare DM to original
        print(f"\nComparing DM to original:")
        res = wilcoxon(R_dm[:, 1], R_dm[:, 0], alternative='greater', nan_policy='omit')
        print(f"{list(SC_dm.keys())[1]} vs {list(SC_dm.keys())[0]}: statistic={res.statistic:.2f}, p-value={res.pvalue:.0e}")
        
        res = wilcoxon(R_dm[:, 3], R_dm[:, 2], alternative='greater', nan_policy='omit')
        print(f"{list(SC_dm.keys())[3]} vs {list(SC_dm.keys())[2]}: statistic={res.statistic:.2f}, p-value={res.pvalue:.0e}")
        
        # Count significant correlations for DM
        sig_counts_dm = (p_dm < 0.05).sum(axis=0)
        print(f"\n=== Number of Significant Age Correlations (p < 0.05) for {exp_id} ===")
        for i, key in enumerate(SC_dm.keys()):
            print(f"{key}: {sig_counts_dm[i]}")
        
        # Minimum p-values (using nanmin to exclude NaNs)
        min_p_dm = np.nanmin(p_dm, axis=0)
        print(f"\n=== Minimum Bonferroni-Corrected p-values for {exp_id} ===")
        for i, key in enumerate(SC_dm.keys()):
            if np.isnan(min_p_dm[i]):
                print(f"{key}: No valid p-values (all NaN)")
            else:
                print(f"{key}: {min_p_dm[i]:.0e}")



# ======= Dataset-specific analyses =======
# Analyze correlations within each individual dataset (before harmonization)
#datasets = ['OASIS-3', 'ADNI2', 'PREVENT-AD']  # Updated PREVENT-AD name

datasets = ['OASIS-3','ADNI2','PREVENT-AD']
print("\n\n===== Within-Dataset Age Correlations (Before Harmonization) =====")
for dataset in datasets:
    for c in ['C','aC']:
        r = np.zeros(aC.shape[1])
        p = np.zeros(aC.shape[1])
        
        age_index = data.loc[(~data.Age.isna()) & (data['SITE']==dataset)].index
        age_values = data.loc[age_index, 'Age'].values
        fea = norm_fea(C if c == 'C' else aC)[age_index, :]

        #fea = (C if c == 'C' else aC)[age_index, :]
        
        for i in range(fea.shape[1]):
            r[i], p[i] = scipy.stats.pearsonr(age_values, fea[:, i])
        
        # Replace NaNs if needed
        if np.isnan(p).any():
            print(f"Warning: Found {np.isnan(p).sum()} NaN p-values in {dataset} {c}. Replacing with 1.0")
            p = np.nan_to_num(p, nan=1.0)
        
        if np.isnan(r).any():
            print(f"Warning: Found {np.isnan(r).sum()} NaN r-values in {dataset} {c}. Replacing with 0.0")
            r = np.nan_to_num(r, nan=0.0)
        
        R = np.abs(r)
        mean_r = np.nanmean(R, axis=0)
        std_r = np.nanstd(R, axis=0)
        
        p_corrected = p * 3570  # Bonferroni correction
        
        print(f"\n{c} with {dataset} in age")
        print(f"Mean |r|: {mean_r:.2f} (±{std_r:.2f})")
        print(f"Min p-value: {np.nanmin(p_corrected):.0e}")
        print(f"Significant connections: {(p_corrected < 0.05).sum()}")



def print_table_format_results(data, C_type, base_SC, results, datasets):
    """
    Print results in a format that makes it easy to fill in the table.
    C_type should be 'C' or 'aC' to match the table rows.
    """
    print(f"\n\n===== TABLE FORMAT OUTPUT FOR {C_type} =====")
    
    # First section: Individual datasets before harmonization
    print("\n----- INDIVIDUAL DATASETS -----")
    for dataset in datasets:
        # Calculate dataset-specific stats
        age_index = data.loc[(~data.Age.isna()) & (data['SITE']==dataset)].index
        if len(age_index) == 0:
            print(f"{dataset}: No valid data")
            continue
            
        age_values = data.loc[age_index, 'Age'].values
        fea = norm_fea(C if C_type == 'C' else aC)[age_index, :]
        
        r = np.zeros(fea.shape[1])
        p = np.zeros(fea.shape[1])
        
        for i in range(fea.shape[1]):
            r[i], p[i] = scipy.stats.pearsonr(age_values, fea[:, i])
        
        # Handle NaNs
        r = np.nan_to_num(r, nan=0.0)
        p = np.nan_to_num(p, nan=1.0)
        
        R_abs = np.abs(r)
        mean_r = np.nanmean(R_abs)
        std_r = np.nanstd(R_abs)
        
        p_corrected = p * 3570
        sig_count = (p_corrected < 0.05).sum()
        min_p = np.nanmin(p_corrected)
        
        print(f"{dataset}:")
        print(f"  |r| mean ± std: {mean_r:.2f} ± {std_r:.2f}")
        print(f"  # p < 0.05: {sig_count}")
        print(f"  min p: {min_p:.0e}")
    
    # Second section: Combined dataset with no harmonization
    print("\n----- COMBINED DATASET - NO HARMONIZATION -----")
    idx = 0 if C_type == 'C' else 7  # Index for C or aC in base_SC
    key = list(base_SC.keys())[idx]
    
    r_data = base_r[:, idx]
    p_data = base_p[:, idx]
    
    # Handle NaNs
    r_data = np.nan_to_num(r_data, nan=0.0)
    p_data = np.nan_to_num(p_data, nan=1.0)
    
    R_abs = np.abs(r_data)
    mean_r = np.nanmean(R_abs)
    std_r = np.nanstd(R_abs)
    
    sig_count = (p_data < 0.05).sum()
    min_p = np.nanmin(p_data)
    
    print(f"No harmonization ({key}):")
    print(f"  |r| mean ± std: {mean_r:.2f} ± {std_r:.2f}")
    print(f"  # p < 0.05: {sig_count}")
    print(f"  min p: {min_p:.0e}")
    
    # Third section: Distribution Matching
    print("\n----- DISTRIBUTION MATCHING -----")
    for ref_site in ['OASIS', 'PREVENTAD']:
        for use_gender in [True, False]:
            exp_id = f"DM_{ref_site}_{'withGender' if use_gender else 'noGender'}"
            if exp_id not in results:
                continue
                
            SC_keys = results[exp_id]['SC_keys']
            idx = 1 if C_type == 'C' else 3  # Index for DM_C or DM_aC
            
            r_data = results[exp_id]['r'][:, idx]
            p_data = results[exp_id]['p'][:, idx]
            
            # Handle NaNs
            r_data = np.nan_to_num(r_data, nan=0.0)
            p_data = np.nan_to_num(p_data, nan=1.0)
            
            R_abs = np.abs(r_data)
            mean_r = np.nanmean(R_abs)
            std_r = np.nanstd(R_abs)
            
            sig_count = (p_data < 0.05).sum()
            min_p = np.nanmin(p_data)
            
            # Get Wilcoxon test results
            base_idx = 0 if C_type == 'C' else 2  # Index for base C or aC in results
            res = wilcoxon(np.abs(results[exp_id]['r'][:, idx]), 
                         np.abs(results[exp_id]['r'][:, base_idx]), 
                         alternative='greater', nan_policy='omit')
            wilcoxon_p = res.pvalue
            
            gender_text = "w/ sex" if use_gender else "w/o sex"
            print(f"Reference: {ref_site}, Dist. {gender_text}:")
            print(f"  |r| mean ± std: {mean_r:.2f} ± {std_r:.2f}")
            print(f"  Wilcoxon p: {wilcoxon_p:.0e}")
            print(f"  # p < 0.05: {sig_count}")
            print(f"  min p: {min_p:.0e}")
    
    # Fourth section: ComBat
    print("\n----- COMBAT -----")
    combat_indices = [1, 2, 3] if C_type == 'C' else [8, 9, 10]
    combat_labels = ["w/ age, w/ sex", "w/o age, w/ sex", "w/o age, w/o sex"]
    
    for i, (idx, label) in enumerate(zip(combat_indices, combat_labels)):
        r_data = base_r[:, idx]
        p_data = base_p[:, idx]
        
        # Handle NaNs
        r_data = np.nan_to_num(r_data, nan=0.0)
        p_data = np.nan_to_num(p_data, nan=1.0)
        
        R_abs = np.abs(r_data)
        mean_r = np.nanmean(R_abs)
        std_r = np.nanstd(R_abs)
        
        sig_count = (p_data < 0.05).sum()
        min_p = np.nanmin(p_data)
        
        # Get Wilcoxon test results
        base_idx = 0 if C_type == 'C' else 7
        res = wilcoxon(np.abs(base_r[:, idx]), 
                     np.abs(base_r[:, base_idx]), 
                     alternative='greater', nan_policy='omit')
        wilcoxon_p = res.pvalue
        
        print(f"ComBat {label}:")
        print(f"  |r| mean ± std: {mean_r:.2f} ± {std_r:.2f}")
        print(f"  Wilcoxon p: {wilcoxon_p:.0e}")
        print(f"  # p < 0.05: {sig_count}")
        print(f"  min p: {min_p:.0e}")
    
    # Fifth section: CovBat
    print("\n----- COVBAT -----")
    covbat_indices = [4, 5, 6] if C_type == 'C' else [11, 12, 13]
    covbat_labels = ["w/ age, w/ sex", "w/o age, w/ sex", "w/o age, w/o sex"]
    
    for i, (idx, label) in enumerate(zip(covbat_indices, covbat_labels)):
        r_data = base_r[:, idx]
        p_data = base_p[:, idx]
        
        # Handle NaNs
        r_data = np.nan_to_num(r_data, nan=0.0)
        p_data = np.nan_to_num(p_data, nan=1.0)
        
        R_abs = np.abs(r_data)
        mean_r = np.nanmean(R_abs)
        std_r = np.nanstd(R_abs)
        
        sig_count = (p_data < 0.05).sum()
        min_p = np.nanmin(p_data)
        
        # Get Wilcoxon test results
        base_idx = 0 if C_type == 'C' else 7
        res = wilcoxon(np.abs(base_r[:, idx]), 
                     np.abs(base_r[:, base_idx]), 
                     alternative='greater', nan_policy='omit')
        wilcoxon_p = res.pvalue
        
        print(f"CovBat {label}:")
        print(f"  |r| mean ± std: {mean_r:.2f} ± {std_r:.2f}")
        print(f"  Wilcoxon p: {wilcoxon_p:.0e}")
        print(f"  # p < 0.05: {sig_count}")
        print(f"  min p: {min_p:.0e}")





# Add this after all calculations but before saving results
datasets = ['OASIS-3', 'ADNI2', 'PREVENT-AD']

# Print table format for C (structural connectivity)
print_table_format_results(data, 'C', base_SC, results, datasets)

# Print table format for aC (adjusted connectivity)
print_table_format_results(data, 'aC', base_SC, results, datasets)




# Save results to file
# Check for NaNs in base_r and base_p and replace if needed
if np.isnan(base_p).any():
    print(f"Warning: Found {np.isnan(base_p).sum()} NaN p-values in base_p. Replacing with 1.0")
    base_p = np.nan_to_num(base_p, nan=1.0)

if np.isnan(base_r).any():
    print(f"Warning: Found {np.isnan(base_r).sum()} NaN r-values in base_r. Replacing with 0.0")
    base_r = np.nan_to_num(base_r, nan=0.0)

# Combine base results with DM results
base_results = {
    'base_harmonization': {
        'r': base_r,
        'p': base_p,
        'SC_keys': list(base_SC.keys())
    }
}
all_results = {**base_results, **results}

import pickle
with open(f'./Results_revision/age_correlation_results.pkl', 'wb') as f:
    pickle.dump(all_results, f)

# Also save a summary DataFrame for easier analysis
summary_rows = []

# Add base harmonization results
R_base = np.abs(base_r)
mean_r_base = np.nanmean(R_base, axis=0)
std_r_base = np.nanstd(R_base, axis=0)
sig_counts_base = (base_p < 0.05).sum(axis=0)

for i, key in enumerate(base_SC.keys()):
    summary_rows.append({
        'experiment': 'base_harmonization',
        'matrix': key,
        'mean_abs_r': mean_r_base[i],
        'std_abs_r': std_r_base[i],
        'sig_connections': sig_counts_base[i],
        'min_p_value': np.nanmin(base_p[:, i]) if not np.all(np.isnan(base_p[:, i])) else float('nan')
    })

# Add DM results
for exp_id in results:
    R_dm = np.abs(results[exp_id]['r'])
    mean_r_dm = np.nanmean(R_dm, axis=0)
    std_r_dm = np.nanstd(R_dm, axis=0)
    sig_counts_dm = (results[exp_id]['p'] < 0.05).sum(axis=0)
    
    for i, key in enumerate(results[exp_id]['SC_keys']):
        min_p_val = np.nanmin(results[exp_id]['p'][:, i]) if not np.all(np.isnan(results[exp_id]['p'][:, i])) else float('nan')
        summary_rows.append({
            'experiment': exp_id,
            'matrix': key,
            'mean_abs_r': mean_r_dm[i],
            'std_abs_r': std_r_dm[i],
            'sig_connections': sig_counts_dm[i],
            'min_p_value': min_p_val
        })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv('./Results_revision/age_correlation_summary.csv', index=False)

# Also save the absolute correlation values for potential visualization
abs_corr_results = {}
# Save base harmonization absolute correlations
abs_corr_results['base_harmonization'] = R_base

# Save DM absolute correlations
for exp_id in results:
    abs_corr_results[exp_id] = np.abs(results[exp_id]['r'])

with open(f'./Results_revision/abs_corr_values.pkl', 'wb') as f:
    pickle.dump(abs_corr_results, f)

print("\nAnalysis complete! Results saved to ./Results_revision/ directory.")




#correlation_results, summary = correlate_connectivity_with_age(conn_mat, dm_har_fea, combat_har_fea, covbat_har_fea, data)