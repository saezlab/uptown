import pandas as pd
import numpy as np
from scipy import optimize, interpolate
import os
import matplotlib.pyplot as plt


def find_files(root_dir, filename):
    matches = []
    for dirpath, dirnames, files in os.walk(root_dir):
        for file in files:
            if file == filename:
                matches.append(os.path.join(dirpath, file))
    return matches

# Uso de la función
root_dir = './decryptm'  # reemplaza esto con tu directorio raíz
filename = 'ratio_single-site_MD.tsv'
matches = find_files(root_dir, filename)

experimental_summary = pd.ExcelFile('./decryptm_old/Experiment_summary.xlsx')

# get PTMs tab
ptm_summary = experimental_summary.parse('PTMs')
ptm_summary.columns
# get columns Conc0-Conc10 and unique values
concentrations = ptm_summary.iloc[:, 19:29][(ptm_summary['Cell line'] == 'A431') & (ptm_summary['Drug'].isin(['Afatinib', 'Dasatinib', 'Gefitinib']))].dropna(axis=0)
concentrations = concentrations.values.flatten()
# get non-numeric values
concentrations = [x for x in concentrations if isinstance(x, (int, float))]

concentrations = np.unique(concentrations)
# arrange them in ascending order
concentrations = np.sort(concentrations)

# create dict matching sample-02 to sample-11 to concentrations
concentrations_dict = {}
for i in range(2, 12):
    sample_col = f"sample-{str(i).zfill(2)}"
    concentrations_dict[sample_col] = concentrations[i - 2]

# per path, read file
df_list = []
for match in matches:
    df = pd.read_table(match, sep='\t')

    # get id from splitting path by /
    id = match.split('/')[2]

    columns_of_interest = [f"sample-{str(i).zfill(2)}" for i in range(2, 12)]

    # remove sample-01 column
    df['runid'] = id
    # rename conc columns
    # remove rows with at least 1 nan in the sample columns
    df.dropna(subset=columns_of_interest, inplace=True)

    df_list.append(df)

# concat all dfs
df = pd.concat(df_list)

df.dropna(axis=1, how='all', inplace=True)

# filter collectri TFs
collectri_tfs = np.unique(pd.read_table('collectri.tsv', sep = '\t').source.values)



# filter df by collectri tfs

# separate runid column in several columns by a separator
df[['analysistype', 'cell_line', 'drug', 'time', 'rep']] = df['runid'].str.split('_', expand=True)


df.columns
# ['Spectrum', 'Spectrum File', 'Peptide',
#        'Extended Peptide', 'Prev AA', 'Next AA', 'Peptide Length', 'Charge',
#        'Retention', 'Observed Mass', 'Calibrated Observed Mass',
#        'Observed M/Z', 'Calibrated Observed M/Z', 'Calculated Peptide Mass',
#        'Calculated M/Z', 'Delta Mass', 'Expectation', 'Hyperscore',
#        'Nextscore', 'PeptideProphet Probability',
#        'Number of Enzymatic Termini', 'Number of Missed Cleavages',
#        'Protein Start', 'Protein End', 'Intensity', 'Assigned Modifications',
#        'Observed Modifications', 'M:15.9949', 'M:15.9949 Best Localization',
#        'STY:79.9663', 'STY:79.9663 Best Localization', 'Purity', 'Is Unique',
#        'Protein', 'Entry Name', 'Protein Description',
#        'Mapped Genes', 'Mapped Proteins', 'Quan Usage',
#        'Modified Peptide', 'Protein ID', 'Gene', 'runid', 'analysistype', 'cell_line', 'drug', 'time', 'rep']

df.drop(columns=['ProteinID', 'Peptide', 'SequenceWindow', 'MaxPepProb',
       'ReferenceIntensity'], inplace=True)

df_columns = [
    #    'Modified Peptide',
    #    'Retention', # same proteins can elute at different times, that's why we keep this column, to prevent duplicates
    #    'Assigned Modifications',
       'Index', 'Gene', 'runid',
       'analysistype', 'cell_line', 'drug', 'time', 'rep']

df_columns = ['Index', 'Gene', 'runid', 'analysistype', 'cell_line', 'drug',
       'time', 'rep']

# pass to long format
longf_df = pd.melt(df, id_vars=df_columns,
             var_name='sample', value_name='ratio')

longf_df['conc'] = longf_df['sample'].map(concentrations_dict)

longf_df.to_csv('detected_peptides_long_MD.csv', index=False)

phospho_prots = longf_df[['drug', 'Gene']].drop_duplicates()

phospho_prots.to_csv('phospho_prots.tsv', sep='\t', index=False)


longdf_martin = pd.read_csv('detected_peptides_long.csv')

# per Modified Peptide, Gene, drug, cell_line and rep, get the number of rows in the longdf_martin
longdf_size = longf_df.groupby(df_columns).size().reset_index(name='counts')
# get entries with counts distinct to 10
longdf_size[longdf_size['counts'] != 10]

# count nmber of rows per peptide, gene, drug, cell_line and rep in the df object

longf_df = longf_df[longf_df['Gene'].isin(collectri_tfs)]

# plot ratios distribution
plt.figure()
plt.hist(longf_df['ratio'], bins=100)
plt.show()

# from https://github.com/kusterlab/decryptM

def logistic_model(x, ec50, slope, top, bottom):
    """
    Logistic model to fit the drug response data.

    Parameters
    ----------
    x: array-like
        Drug concentrations in log space
    ec50: float
        Inflection point in log space of the sigmoid value
    slope: float
        slope of the transition between top and bottom
    top: float
        top asymptote
    bottom: float
        bottom asymptote

    Returns
    -------
    y: array-like
        Drug response
    """
    # predict y with given parameters using the standard model
    return (top - bottom) / (1 + 10 ** (slope * (x - ec50))) + bottom

def fit_logistic_function(y_data, x_data, x_interpolated=None, curve_guess=None, curve_bounds=None, max_opt=None):
    """
    fit_logistic_function(y_data, x_data, x_interpolated=None, max_opt=None)

    This function is fitting a 4 parameter logistic function to the given y & x data.

    Parameters
    ----------
    y_data : pd.Series
        A Series object, containing the observed DMSO normalized ratios (including the DMSO ratio)
    x_data : array_like
        An array, containing the used drug concentrations in Log space
    x_interpolated : array_like, optional
        If None (default), the fit will be performed on the pure y&x data. If object is an array,
        containing drug concentrations, then the logistic fit is performed on those x_interpolated
        values. To obtain y_interpolated values that are corresponding to the x_interpolated, a
        linear interpolation is performed. The values of x_interpolated should be within within the
        x_data range to not extrapolate the data.
    curve_guess : array_like, optional
        If None (default), the fit will use no initial parameter guesses.
    curve_bounds : tuple of array_like, optional
        If None (default), the fit will use no bounds for the parameter values.
        Else provide the following structure ([min1, ..], [max1, ..]).
    max_opt : int, optional
        If None (default), the fit will use scipy default estimation of a maximum number of fit
        iterations before the fit is terminated. Otherwise, provide a integer to specify the maximum
        number yourself to save time or increase number of successful fits in longer time.

    Returns
    -------
    out : tuple
        An tuple object containing the fitted parameter, R2, and the parameter standard errors in
        the following order: ('Log EC50', 'Curve slope', 'Curve top', 'Curve bottom', 'R2', 'Curve RMSE',
        'Log EC50 error', 'Curve slope error', 'Curve top error', 'Curve bottom error')
    """
    # fast return
    if any(np.isnan(y_data)):
        return 10 * (np.nan,)

    # Fit logistic model
    try:
        y_data = y_data.values

        if x_interpolated is not None:
            f = interpolate.interp1d(x_data, y_data, kind='linear')
            popt, pcov = optimize.curve_fit(logistic_model, x_interpolated, f(x_interpolated),
                                            p0=curve_guess, bounds=curve_bounds, max_nfev=max_opt)
            perr = np.sqrt(np.diag(pcov))

        else:
            popt, pcov = optimize.curve_fit(logistic_model, x_data, y_data, p0=curve_guess, maxfev=5000)
            perr = np.sqrt(np.diag(pcov))

        # Calculate R2 and RMSE
        ss_res = np.sum((y_data - logistic_model(x_data, *popt)) ** 2)
        ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
        r2 = 1 - ((ss_res + 1e-10) / (ss_tot + 1e-10))
        rmse = np.sqrt(ss_res/len(y_data))

        return (*popt, r2 if r2 >= 0 else 0, rmse, *perr)

    except RuntimeError:
        return 10 * (np.nan,)

# using the concentrations dict, create a new column with the concentration for each sample


longf_df
# per gene, drug and cell_line, fit logistic curve
# create empty df
fit_params = pd.DataFrame(columns=['Index', 'Gene', 'Drug', 'Cell_line', 'rep', 'log_ec50', 'slope', 'top', 'bottom', 'rsq', 'rmse', 'log_ec50_error', 'slope_error', 'top_error', 'bottom_error'])


# get unique combinations of gene, drug and cell_line
unique_combinations = longf_df[['Index', 'Gene', 'drug', 'cell_line', 'rep']].drop_duplicates()


for comb in unique_combinations.values:
    psite, gene, drug, cell_line, rep = comb
    x_data = np.log10(longf_df[(longf_df['Index'] == psite) & (longf_df['Gene'] == gene) & (longf_df['drug'] == drug) & (longf_df['cell_line'] == cell_line) & (longf_df['rep'] == rep)]['conc'].values)
    y_data = longf_df[(longf_df['Index'] == psite) & (longf_df['Gene'] == gene) & (longf_df['drug'] == drug) & (longf_df['cell_line'] == cell_line) & (longf_df['rep'] == rep)]['ratio']
    log_ec50, slope, top, bottom, r2, rmse, log_ec50_error, slope_error, top_error, bottom_error = fit_logistic_function(x_data=x_data, y_data=y_data, curve_guess=[np.median(x_data), 1, 0, 0])


    # fit_params = fit_params.append({'Index': psite, 'Gene': gene, 'Drug': drug, 'Cell_line': cell_line, 'rep': rep, 'log_ec50': log_ec50, 'slope': slope, 'top': top, 
    #                                 'bottom': bottom, 'rsq': r2, 'rmse': rmse, 'log_ec50_error': log_ec50_error, 'slope_error': slope_error, 
    #                                 'top_error': top_error, 'bottom_error': bottom_error}, ignore_index=True)
    
    # # plot the fit
    if r2 > 0.8:
        plt.figure()
        plt.scatter(x_data, y_data)
        plt.plot(x_data, logistic_model(x_data, log_ec50, slope, top, bottom))
        plt.title(f"{psite} {drug} {cell_line}")
        # x axis as log
        # plt.savefig(f"{gene}_{drug}_{cell_line}.png")
        plt.show()



# save fitted curves and parameters
fit_params.to_csv('fit_params_peptide_rep_MD.csv', index=False)


# 131 genes with peptides with at least r2 > 0.8
len(fit_params[fit_params['rsq'] > 0.8].unique())



