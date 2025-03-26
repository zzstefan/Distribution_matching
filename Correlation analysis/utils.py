####In this script, we will include some function which will be frequently used;

import numpy as np
import pandas as pd
import scipy
import neuroHarmonize as nh
from scipy.stats import gamma
from covbat_fixed import covbat as covbat_fixed


def read_sc(dataframe, C_measure):
    baseline = dataframe
    sc = np.zeros((len(baseline), int(85 * 84 / 2)))
    if C_measure == 'C':
        for i in range(len(baseline)):
            tmp = scipy.io.loadmat(baseline.loc[i, 'conn_mat'])['C']
            sc[i, :] = tmp[np.triu_indices(85, k=1)]
    else:
        for i in range(len(baseline)):
            tmp = scipy.io.loadmat(baseline.loc[i, 'conn_mat'])['aC']
            sc[i, :] = tmp[np.triu_indices(85, k=1)]

    return sc


def read_DM_C(dataframe, C_measure):
    DM_C = np.zeros((len(dataframe), int(85*84/2)))
    if C_measure == 'C':
        for i in range(len(dataframe)):
            tmp = scipy.io.loadmat(dataframe.loc[i, 'DM_C'])['DM_C']
            DM_C[i, :] = tmp[np.triu_indices(85, k=1)]
    else:
        for i in range(len(dataframe)):
            tmp = scipy.io.loadmat(dataframe.loc[i, 'DM_aC'])['DM_aC']
            DM_C[i, :] = tmp[np.triu_indices(85, k=1)]

    return DM_C


def outliner(dataframe, C_measure):
    sc = read_sc(dataframe, C_measure)
    Mean_sc = np.mean(sc, axis=1)
    Median_sc = np.median(sc, axis=1)
    dataframe['mean_sc'] = Mean_sc.reshape(-1, 1)
    dataframe['median_sc'] = Median_sc.reshape(-1, 1)

    index_drop = np.where((Mean_sc > 40000) | (Mean_sc < 100) | (Median_sc > 30000) | (Median_sc < 100))[0]
    index_keep = np.where(~((Mean_sc > 40000) | (Mean_sc < 100) | (Median_sc > 30000) | (Median_sc < 100)))[0]
    return index_drop, index_keep



def vec_to_sym(vec, n):
    # vec is the vector to be transformed back to the matrix
    # n is the shape of the symmetric matrix
    # N is the number of subjects
    tmp = np.zeros((n, n))
    inds = np.triu_indices_from(tmp, k=1)
    harmonized_matrix = np.zeros((vec.shape[0], n, n))

    for i in range(vec.shape[0]):
        harmonized_matrix[i, inds[0], inds[1]] = vec[i, :]
        harmonized_matrix[i, inds[1], inds[0]] = vec[i, :]

    return harmonized_matrix


def harmonize(features: object, dataframe: object, use_GAMs: object) -> object:
    dataframe['SEX'] = (dataframe['Sex'] == 'M').astype(int)
    covars = dataframe[['SITE', 'Age', 'SEX']]
    if use_GAMs == 1:
        model, harmonized_features = nh.harmonizationLearn(features, covars,
                                                           smooth_terms=['Age'],
                                                           smooth_term_bounds=(
                                                               np.floor(np.min(dataframe.Age)),
                                                               np.ceil(np.max(dataframe.Age))))
    else:
        model, harmonized_features = nh.harmonizationLearn(features, covars)
    return model, harmonized_features


def harmonize_noAge(features: object, dataframe: object) -> object:
    dataframe['SEX'] = (dataframe['Sex'] == 'M').astype(int)
    covars = dataframe[['SITE', 'SEX']]
    model, harmonized_features = nh.harmonizationLearn(features, covars)
    return model, harmonized_features




def harmonize_noAge_noSex(features: object, dataframe: object) -> object:
    covars = dataframe[['SITE']]
    model, harmonized_features = nh.harmonizationLearn(features, covars)
    return model, harmonized_features


# def apply_covbat2(data, batch, age=None, sex=None):
    """
    Apply CovBat with flexible covariate options: age and sex, sex only, or no covariates.
    
    Parameters:
    - data: np.ndarray, shape (n_subjects, n_features), imaging features
    - batch: np.ndarray, shape (n_subjects,), site/scanner IDs
    - age: np.ndarray, shape (n_subjects,), continuous covariate (optional, default=None)
    - sex: np.ndarray, shape (n_subjects,), categorical covariate (optional, default=None)
    
    Returns:
    - harmonized_data: np.ndarray, shape (n_features, n_subjects), harmonized data
    """
    # Transpose data to subjects x features for covbat
    
    # Determine covariate setup based on inputs
    data_T = data.T
    if age is not None and sex is not None:
        # Case 1: Both age and sex provided
        covariates = pd.DataFrame({
            'age': age,
            'sex': sex
        })
    elif sex is not None:
        # Case 2: Only sex provided
        covariates = pd.DataFrame({
            'sex': sex
        })
    else:
        # Case 3: No covariates (age and sex are None)
        covariates = None
    
    # Apply CovBat
    harmonized_data_T = covbat(
        data=data_T,
        batch=batch,
        numerical_covariates=covariates  # Will be None if no covariates
    )
    
    return harmonized_data_T.T



# def apply_covbat(data, batch, age=None, sex=None):
    """
    Apply CovBat with flexible covariate options: age and sex, sex only, or no covariates.
    
    Parameters:
    - data: np.ndarray, shape (n_subjects, n_features), imaging features
    - batch: np.ndarray or list, shape (n_subjects,), site/scanner IDs
    - age: np.ndarray, shape (n_subjects,), continuous covariate (optional, default=None)
    - sex: np.ndarray, shape (n_subjects,), categorical covariate (optional, default=None)
    
    Returns:
    - harmonized_data: np.ndarray, shape (n_features, n_subjects), harmonized data
    """
    from covbat_fixed import covbat as covbat_fixed
    
    # Ensure all data is in the right format
    data_T = data.T
    
    # Create proper model DataFrame
    model = pd.DataFrame({'batch': batch})
    
    # Add covariates to model if provided
    if age is not None:
        model['Age'] = age.astype(float)  # Ensure age is numeric
    
    if sex is not None:
        # Convert sex to numeric (0/1) if it's not already
        if isinstance(sex[0], str):
            # Assuming sex is either "M"/"F" or "male"/"female"
            sex_numeric = np.array([1 if s in ['M', 'male'] else 0 for s in sex])
            model['Sex'] = sex_numeric
        else:
            model['Sex'] = sex
    
    # Define numerical covariates list
    numerical_covariates = []
    if age is not None:
        numerical_covariates.append(1)
    
    # Apply CovBat with properly formatted model and covariates
    harmonized_data_T = covbat_fixed(
        data=data_T,
        batch=batch,
        model=model,
        numerical_covariates=numerical_covariates
    )
    
    return harmonized_data_T.T



def norm_fea(fea):
    new_fea = np.zeros((fea.shape))
    mean_fea = np.mean(fea, axis=0)
    for i in range(fea.shape[1]):
        if mean_fea[i] == 0:
            new_fea[:, i] = 0
        else:

            new_fea[:, i] = 1 - np.exp(-fea[:, i] / mean_fea[i])
    return new_fea


def recon_full(middle_result, new_data, nonZero_col):
    full_data = np.zeros((len(new_data), 3570))
    full_data[:, nonZero_col] = middle_result
    return full_data

### weight equals (# sample points where c=0)/(# all sample points)
# def mixed_distribution(x, amplitude, shape, loc, scale):
#     return amplitude * np.where(x == 0, 1, 0) + (1 - amplitude) * gamma.pdf(x, shape, loc, scale=scale)


