### All questions, comments, or concerns about the following code can be directed to michelle.mastrianni@fda.hhs.gov.
#
### This notebook takes an input dataset which contains:
#
# (1) A column containing true cancer labels per image ('label').
# (2) A column containing AI scores for each image ('AIscores'). These can be on any scale.
# (3) *at least one of* (a) A column containing radiologist BIRADS scores corresponding to the screening image ('BIRADS') -- 0, 1, or 2, with 0 meaning 'recall' and 1 and 2 'negative', AND/OR
#                       (b) A column containing radiologist scores on some finer scale, e.g., from reader studies, ('radscores').
# (4) A column containing reader IDs per radiologist ('radID') (OPTIONAL)
# (5) A column containing image IDs ('imageID') (OPTIONAL)
# (6) A column containing BIRADS density scores (1-4 or A-D). (OPTIONAL)
#
### The purpose of this notebook is to:
#
# (Note: the following applies on a per-reader and/or per-density basis, depending on user input specifications.)
# (1) (a) Calculate the radiologist per-image sensitivity and specificity (if BI-RADS is used), AND/OR
#     (b) Plot the empirical radiologist ROC curve (if finer scale of radiologist scores is used) and calculate its AUC, AND
#     (c) (if applicable) plot empirical radiologist ROC curves for individual radiologists.
# (2) Plot the empirical AI ROC curve and calculate its AUC
# (3) Fit a parametric ROC curve for the AI diagnostic
# (4) (a) using parameters from AI ROC curve, fit a parametric radiologist ROC curve to single (FPF, TPF) radiologist point (if BI-RADS is used), AND/OR
#     (b) fit a parametric ROC curve for the radiologist diagnostic (if finer scale of radiologist scores is used)
# (5) Compare AI scores with radiologist screening assessment to get correlation measures
# (6) (if applicable) get correlation results between individual radiologists.
# (6) Generate CSV file (and some plots) containing empirical rule-out data (for example, calculating the empirical performance of a radiologist under a rule-out, rule-in, or combined scenario, as well as the estimated performance via copulas theory).
#
# The expected file outputs are as follows. Outputs are saved in folders outputs/outputs_ID_i_density_j where i represents the reader ID and j represents the BIRADS density category. If not conditioning on a specified category, i or j is "full".
#
# (1) 'summary.csv'. Saves number of images and number of cancer images.
# (2) 'AI_score_histogram.png'. This plots a score histogram for the AI. This is the only image file that is saved without the data used to generate it, as it is formed from the AI scores themselves, which our research group should not have access to.
# (3) 'empirical_fpr_tpr_AI_All.csv'. This gives the fpr and tpr values for the empirical AI ROC curve at different thresholds.
# (4) (IF 'BIRADS' column): 'empirical_fpr_tpr_birads.csv'. This simply outputs the radiologist performance (FPR, TPR) using BI-RADS scores (0 = yes, 1 or 2 = no)
# (5) (IF 'radscores' column): 'empirical_fpr_tpr_rad_All.csv'. This gives the fpr and tpr values for the empirical radiologist ROC curve at different thresholds.
# (6) 'empirical_ROC_curves_All.png': This plots (3), (4), and (5).
# (7) SEVERAL FILES (IF 'radID' column): 'empirical_fpr_tpr_birads_i.csv' where i ranges from 1 to the total number of individual radiologists. Same idea as (3) but grouped by each individual reader.
# (8) SEVERAL FILES (IF 'radID' column): 'empirical_fpr_tpr_rad_i.csv' where i ranges from 1 to the total number of individual radiologists. Same idea as (4) but grouped by each individual reader.
# (9) SEVERAL FILES (IF 'radID' column): 'empirical_ROC_curves_i.png': plots (3), (7) and (8) for each reader i.
# (10) 'deming_AI_line.csv'. This outputs (x, y) values to plot the Deming regression line for the AI fit.
# (11) 'deming_AI_points.csv'. This outputs the transformed empirical AI (FPR, TPR) points in Deming space.
# (12) 'Deming_AI.png'. This plots (10) and (11).
# (13) 'AI_FPR_TPR_fit.csv'. This gives the fpr and tpr values for the smooth/fitted AI ROC curve at different thresholds.
# (14) 'empirical_fitted_AI.png'. This plots (13).
# (15) (IF 'radscores' column): deming_rad_line.csv. This outputs (x, y) values to plot the Deming regression line for the radiologist (if applicable).
# (16) (IF 'radscores' column): deming_rad_points.csv. This outputs the transformed empirical radiologist (FPR, TPR) points in Deming space.
# (17) (IF 'radscores' column): Deming_rad.png. This plots (15) and (16).
# (18) (IF 'radscores' column): 'rad_FPR_TPR_fit_continuous.csv'. This gives the fpr and tpr values for the smooth/fitted radiologist ROC curve at different thresholds.
# (19) (IF 'radscores' column): 'empirical_fitted_rad.png'. This plots (18).
# (20) (IF 'BIRADS' column): 'rad_FPR_TPR_fit_birads.csv'. This gives the fpr and tpr values for the smooth/fitted radiologist curve to the (FPR, TPR) point from BIRADS scores assuming the same mean-to-sigma ratio as the AI ROC curve.
# (21) (IF 'BIRADS' column): 'smooth_AI_rad.png'. This plots (20).
# (22) 'a_b_values.csv': stores the a and b values necessary for calculating mu and sigma for the various fitted ROC curves.
# (23) (IF 'BIRADS' column): 'correlation_results_BIRADS_All.csv'. This computes Pearson's rho and Kendall's tau correlations between AI scores and radiologist 0-1-2 BIRADS (as well as 0/1 binary scores)
# (24) (IF 'radscores' column): 'correlation_results_continuous_All.csv'. This computes Pearson's rho and Kendall's tau correlations between AI scores and continuous radiologist scores.
# (25) SEVERAL FILES (IF 'BIRADS' column and 'radID' column): 'correlation_results_BIRADS_i.csv' where i ranges from 1 to total number of radiologists. Same as 'correlation_results_BIRADS' but for each individual reader.
# (26) SEVERAL FILES (IF 'radscores' column and 'radID' column): 'correlation_results_continuous_i.csv' where i ranges from 1 to total number of radiologists. Same as 'correlation_results_continuous' but for each individual reader.
# (27) (IF 'BIRADS' column, 'imageID' column, and 'radID' column): 'all_correlation_results_birads.csv': saves pairwise correlation results between pairs of readers, for all readings (BIRADS case)
# (28) (IF 'BIRADS' column, 'imageID' column, and 'radID' column): 'diseased_correlation_results_birads.csv': saves pairwise correlation results between pairs of readers, for diseased cases (BIRADS case)
# (29) (IF 'BIRADS' column, 'imageID' column, and 'radID' column): 'non_diseased_correlation_results_birads.csv': saves pairwise correlation results between pairs of readers, for non-diseased cases (BIRADS case)
# (30) (IF 'radscores' column, 'imageID' column, and 'radID' column): 'all_correlation_results_continuous.csv': saves pairwise correlation results between pairs of readers, for all readings (finer radiologist scores case)
# (31) (IF 'radscores' column, 'imageID' column, and 'radID' column): 'diseased_correlation_results_continuous.csv': saves pairwise correlation results between pairs of readers, for diseased cases (finer radiologist scores case)
# (32) (IF 'radscores' column, 'imageID' column, and 'radID' column): 'non_diseased_correlation_results_continuous.csv': saves pairwise correlation results between pairs of readers, for non-diseased cases (finer radiologist scores case)
# (33) (IF 'BIRADS' column): 'AI_FPF_rad_performance_birads.csv'. This outputs the empirical radiologist performance under rule-out and rule-in at different rule-out/rule-in AI FPF thresholds.
# (34) (IF 'BIRADS' column): several figures ('rule_out_{FPF_ruleout}_rule_in_{FPF_rulein}_birads'.png) plotting results from the csv file. About 5 such figures will be generated as an example.
# (35) (IF 'radscores' column): 'AI_FPF_rad_performance_continuous.csv'. Same as above, but does so for various working thresholds for the radiologist (assuming continuous radiologist scores). Also outputs data assuming that the radiologist changes his/her operating threshold under rule-out slightly (this is the idea behind the "w" values).
# (36) (IF 'radscores' column): several figures ('rule_out_{FPF_ruleout}_rule_in_{FPF_rulein}_continuous'.png) plotting results from the csv file. About 5 such figures will be generated as an example.


### Import packages

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
import math
from scipy.stats import norm, expon, multivariate_normal, kendalltau, pearsonr, beta
from scipy.optimize import minimize, minimize_scalar
from sklearn.metrics import roc_curve, auc
from scipy.odr import Model, RealData, ODR
from scipy.integrate import quad
import os
import sys
from itertools import combinations
import subprocess
import tabulate
import warnings

warnings.filterwarnings('ignore')

################## Read in data

### The data from the git repo (https://github.com/michellemastrianni/rule-out-example) is example data, randomly generated.
### example_data.csv contains only AI scores and radiologist BIRADS in addition to true cancer labels.
### example_data_with_rad_scores.csv contains AI scores, radiologist BIRADS, radiologist scores on a finer scale, radiologist IDs, and true cancer labels.
### example_data_with_rad_and_image.csv contains AI scores, radiologist BIRADS, radiologist scores on a finer scale, radiologist IDs, image IDs (to obtain correlations between radiologists on the same images), and true cancer labels.
### example_data_with_density.csv contains all of the above and a density column.

### When running your own data:
### (1) First, make sure column names appropriately match the ones in the script (see first cell for details).

if len(sys.argv) < 2:
    print("Usage: python script.py")
    sys.exit(1)
filename = sys.argv[1]
df_orig = pd.read_csv(filename)


def full_analysis(df, ID, density):

    ### Print warnings if missing or unexpected data
    def check_for_nan(df, column_name):
        normal_values = df[df[column_name].notnull()]
        nan_values = df[column_name].isnull()
        if nan_values.sum() > 0:
            print(f"Warning: Found {nan_values.sum()} NaN value(s) in column '{column_name}'.")

    def check_for_unexpected_values(df, column_name):
        unexpected_values = df[~df[column_name].apply(lambda x: isinstance(x, (int, float)))][column_name].tolist()
        if unexpected_values:
            print(f"Warning: Found non-numeric value(s) {unexpected_values} in column '{column_name}'.")

    def check_birads(df):
        valid_values = [0, 1, 2, 'A', 'N', 'B']
        invalid_values = df[~df['BIRADS'].isin(valid_values)]['BIRADS'].tolist()
        if invalid_values:
            print(f"Warning: Found invalid BIRADS value(s) {invalid_values}.")

    #Check for insufficient data. Move to next reader/density group if not enough data to conduct analysis.
    non_diseased_filtered = df[(df['label'] == 0) & (df['AIscores'].notnull())]['AIscores'].count()
    diseased_filtered = df[(df['label'] == 1) & (df['AIscores'].notnull())]['AIscores'].count()
    if non_diseased_filtered < 1 or diseased_filtered < 1:
        print("Insufficient data. Moving to next group.")
        return

    check_for_nan(df, 'label')
    check_for_unexpected_values(df, 'label')
    check_for_nan(df, 'AIscores')
    check_for_unexpected_values(df, 'AIscores')

    if 'BIRADS' in df.columns:
        check_for_nan(df, 'BIRADS')
        check_for_unexpected_values(df, 'BIRADS')
        check_birads(df)
    if 'radscores' in df.columns:
        check_for_nan(df, 'radscores')
        check_for_unexpected_values(df, 'radscores')

    ### Brief summary statistics (total number of images and number of images with cancer)

    df = df.dropna()

    if 'imageID' not in df.columns:
        num_images = len(df)
        num_cancer_images = len(df[df['label'] == 1])
    else:
        num_images = len(df['imageID'].unique())
        num_cancer_images = len(df[df['label'] == 1]['imageID'].unique())

    data = {'total_images': [num_images], 'cancer_images': [num_cancer_images]}
    df_info = pd.DataFrame(data)
    df_info_str = df_info.to_markdown(index=False)
    print(df_info_str)

    # Save the DataFrame to a CSV file

    df_info.to_csv(f'outputs/outputs_ID_{ID}_density_{density}/summary.csv', index=False)

    ### Rescale the data

    min_AI_val = df['AIscores'].min()
    max_AI_val = df['AIscores'].max()
    df['AIscores'] = (df['AIscores'] - min_AI_val) / (max_AI_val - min_AI_val)

    if 'radscores' in df.columns:
        min_rad_val = df['radscores'].min()
        max_rad_val = df['radscores'].max()
        df['radscores'] = (df['radscores'] - min_AI_val) / (max_AI_val - min_AI_val)

    ### Plot re-scaled AI score histograms (separating non-diseased and diseased scores)

    non_diseased_filtered = df[df['label'] == 0]['AIscores']
    diseased_filtered = df[df['label'] == 1]['AIscores']

    plt.figure(figsize=(10, 6))
    freq, bins, _ = plt.hist([non_diseased_filtered, diseased_filtered],
                               bins=30, color=['blue', 'red'],
                               label=['Non-Diseased', 'Diseased'],
                               alpha=0.7, density=True)

    plt.xlabel('AI Scores')
    plt.ylabel('Frequency')
    plt.title('AI Score Histogram')
    plt.legend()
    plt.ylim(0, np.max(freq) + 0.1)
    plt.savefig(f'outputs/outputs_ID_{ID}_density_{density}/AI_score_histogram.png')

    ### Plot empirical ROC curve for AI scores and calculate per-image radiologist FPR/TPR (or empirical ROC if finer radiologist scores are available).

    # The purpose of this function (in ROC_get) is to check if ROC_get is performed on the whole dataframe or on a subgroup corresponding to a specific reader. If specific reader, it gets the ID corresponding to that reader.
    def are_all_same(dfgroup, column_name):
        first_entry = dfgroup[column_name].iloc[0]
        return all(dfgroup[column_name] == first_entry)

    # Function to calculate (FPR, TPR) value for radiologist.
    # If using BIRADS score, we add 'rad' and 'rad_ANB' columns.
    # 'rad' is a binary column keeping track of either a "positive" (1) or "negative" (0) reading.
    # 'rad_ANB' includes some separation between the "1" and "2" BI-RADS designations, with "2" being a slightly more suspicious negative than "1".
    # We do not use 'rad_ANB' later in the code but include it for completeness.
    def cal_birads_point(dfgroup):
        birads_mapping = {
            'N': (0, 0),
            'B': (0, 1),
            'A': (1, 2),
            1: (0, 0),
            2: (0, 1),
            0: (1, 2)
        }
        for index, row in dfgroup.iterrows():
            birads = row['BIRADS']
            if isinstance(birads, str) or birads in [0, 1, 2]:
                rad, rad_ANB = birads_mapping.get(birads, (np.nan, np.nan))
                dfgroup.loc[index, 'rad'] = rad
                dfgroup.loc[index, 'rad_ANB'] = rad_ANB
            else:
                dfgroup.loc[index, 'rad'] = np.nan
                dfgroup.loc[index, 'rad_ANB'] = np.nan

        # Calculate radiologist FPR and TPR point
        TN = len(dfgroup[(dfgroup['rad'] == 0) & (dfgroup['label'] == 0)])
        TP = len(dfgroup[(dfgroup['rad'] == 1) & (dfgroup['label'] == 1)])
        FN = len(dfgroup[(dfgroup['rad'] == 0) & (dfgroup['label'] == 1)])
        FP = len(dfgroup[(dfgroup['rad'] == 1) & (dfgroup['label'] == 0)])
        if TN + FP != 0:
            empirical_fpr_rad_birads = FP / (TN + FP)
        else:
            empirical_fpr_rad_birads = 0
        if TP + FN != 0:
            empirical_tpr_rad_birads = TP / (TP + FN)
        else:
            empirical_tpr_rad_birads = 0
        return empirical_fpr_rad_birads, empirical_tpr_rad_birads

    def ROC_get(dfgroup):
        # Get group name (either 'All' if all readers, or reader ID if a specific reader)
        if 'radID' in dfgroup.columns and condition_reader == 'False':
            if are_all_same(dfgroup, 'radID'):
                groupName = dfgroup['radID'].iloc[0]
            else:
                groupName = 'All'
        else:
            groupName = 'All'

        dfgroup = dfgroup.dropna()
        
        num_cases = len(dfgroup)  # number of cases in the reader group
        empirical_fpr_AI, empirical_tpr_AI, thresholds_binary = roc_curve(dfgroup['label'], dfgroup['AIscores'])
        empirical_auc_AI = auc(empirical_fpr_AI, empirical_tpr_AI)

        if 'BIRADS' in dfgroup.columns:
            empirical_fpr_rad_birads, empirical_tpr_rad_birads = cal_birads_point(dfgroup)

        if 'radscores' in dfgroup.columns:
            empirical_fpr_rad, empirical_tpr_rad, thresholds_binary = roc_curve(dfgroup['label'], dfgroup['radscores'])
            empirical_auc_rad = auc(empirical_fpr_rad, empirical_tpr_rad)

        plt.figure(figsize=(6, 6))

        # Plot AI ROC curve and radiologist FPR/TPR point, and print AUC of ROC curve
        if len(df) == len(dfgroup):  # If we are in the case of 'all' readers, plot the AI ROC along with the radiologist ROC.
            plt.plot(empirical_fpr_AI, empirical_tpr_AI, color='red', label='AI ROC', lw=2)
        if 'BIRADS' in dfgroup.columns:  # Plot single (FPR, TPR) point based on BIRADS
            plt.scatter(empirical_fpr_rad_birads, empirical_tpr_rad_birads, color='green', label='Radiologist BIRADS')
        if 'radscores' in dfgroup.columns:  # Plot empirical ROC curve based on finer radiologist scores
            plt.plot(empirical_fpr_rad, empirical_tpr_rad, color='green', label='Radiologist from Scores')

        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.0])
        plt.xlabel('FPR')
        plt.ylabel('TPR')
        plt.title(f'AI ROC Curve and Radiologist Performance, Reader {groupName}, Cases: {num_cases}')
        plt.legend(loc='lower right')
        plt.grid(True)
        print(f'AI AUC, reader {groupName}: ', empirical_auc_AI)
        if 'radscores' in dfgroup.columns:
            print(f'rad AUC, reader {groupName}: ', empirical_auc_rad)
        if 'BIRADS' in dfgroup.columns:
            print(f'Radiologist FPR, Reader {groupName}: ', empirical_fpr_rad_birads)
            print(f'Radiologist TPR, Reader {groupName}: ', empirical_tpr_rad_birads)

        ### Save files
        plt.savefig(f'outputs/outputs_ID_{ID}_density_{density}/empirical_ROC_curves_{groupName}.png')
        # plt.show()

        empirical_fpr_tpr_AI = pd.DataFrame({
            'empirical_fpr': empirical_fpr_AI,
            'empirical_tpr': empirical_tpr_AI
        })
        empirical_fpr_tpr_AI.to_csv(f'outputs/outputs_ID_{ID}_density_{density}/empirical_fpr_tpr_AI_{groupName}.csv', index=False)

        if 'radscores' in dfgroup.columns:
            empirical_fpr_tpr_rad = pd.DataFrame({
                'empirical_fpr': empirical_fpr_rad,
                'empirical_tpr': empirical_tpr_rad
            })
            empirical_fpr_tpr_rad.to_csv(f'outputs/outputs_ID_{ID}_density_{density}/empirical_fpr_tpr_rad_{groupName}.csv', index=False)

        if 'BIRADS' in dfgroup.columns:
            empirical_fpr_tpr_birads_df = pd.DataFrame({
                'empirical_fpr': [empirical_fpr_rad_birads],
                'empirical_tpr': [empirical_tpr_rad_birads]
            })
            empirical_fpr_tpr_birads_df.to_csv(f'outputs/outputs_ID_{ID}_density_{density}/empirical_fpr_tpr_birads_{groupName}.csv', index=False)
            return empirical_fpr_rad_birads, empirical_tpr_rad_birads

    ROC_get(df)  # Run function on the whole dataframe (all readers)

    ### If individual reader IDs are available, generate empirical FPR/TPR point (if BIRADS) and/or empirical ROC (if finer radiologist scores) for each reader.
    if 'radID' in df.columns:
        random_IDs = df.groupby('radID').ngroup().add(1)  # Replace radiologist IDs with random numbers
        df['radID'] = random_IDs
        grouped = df.groupby('radID')
        grouped.apply(ROC_get)

    #### Transform empirical AI ROC curve using Deming regression

    # Redefine empirical fpr/tpr arrays for AI diagnostic
    empirical_fpr_AI, empirical_tpr_AI, thresholds_AI = roc_curve(df['label'], df['AIscores'])

    def inverse_normal_cdf(p):
        return norm.ppf(p)

    def inverse_exponential_cdf(p):
        return expon.ppf(p)

    def transform_roc_values(empirical_values, cdf_type):
        inverse_funcs = {'normal': inverse_normal_cdf, 'exponential': inverse_exponential_cdf}
        if cdf_type not in inverse_funcs:
            raise ValueError("Invalid 'cdf_type'. Choose 'normal' or 'exponential'.")
        inverse_func = inverse_funcs[cdf_type]
        transformed_values = inverse_func(empirical_values)
        return transformed_values

    def deming_model(B, x):
        return B[0] + B[1] * x

    def fit_deming_regression(x, y):
        model = Model(deming_model)
        data = RealData(x, y)
        odr = ODR(data, model, beta0=[1.0, 1.0])
        output = odr.run()
        intercept, slope = output.beta
        return intercept, slope

    max_first_index = max(np.argmax(empirical_fpr_AI > 0), np.argmax(empirical_tpr_AI > 0))
    min_last_index = min(np.argmax(empirical_fpr_AI >= 1), np.argmax(empirical_tpr_AI >= 1))
    if max_first_index >= min_last_index:
        min_last_index = max_first_index + 1

    # Transform empirical ROC values using inverse normal functions
    transformed_fpr_normal_AI = transform_roc_values(empirical_fpr_AI[max_first_index:min_last_index], 'normal')
    transformed_tpr_normal_AI = transform_roc_values(empirical_tpr_AI[max_first_index:min_last_index], 'normal')

    # Fit Deming regression
    a_normal_AI, b_normal_AI = fit_deming_regression(transformed_fpr_normal_AI, transformed_tpr_normal_AI)

    # Generate x values for Deming regression line
    min_transformed_fpr_AI = np.min(transformed_fpr_normal_AI)
    max_transformed_fpr_AI = np.max(transformed_fpr_normal_AI)
    x_values_normal_AI = np.linspace(min_transformed_fpr_AI - 0.1, max_transformed_fpr_AI + 0.1, 400)
    y_values_normal_AI = deming_model([a_normal_AI, b_normal_AI], x_values_normal_AI)

    # Plot Deming regression line
    plt.figure(figsize=(6, 6))
    plt.plot(x_values_normal_AI, y_values_normal_AI, color='blue', label='Deming Regression Line (Normal)')
    plt.scatter(transformed_fpr_normal_AI, transformed_tpr_normal_AI, color='blue')
    plt.xlabel('Transformed FPR')
    plt.ylabel('Transformed TPR')
    plt.title('Deming Regression for AI - Normal')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'outputs/outputs_ID_{ID}_density_{density}/Deming_AI.png')
    # plt.show()

    # Save x,y values for fitted line and points
    deming_csv_AI_line = pd.DataFrame({
        'x_values_normal': x_values_normal_AI,
        'y_values_normal': y_values_normal_AI
    })
    deming_csv_AI_line.to_csv(f'outputs/outputs_ID_{ID}_density_{density}/deming_AI_line.csv', index=False)

    deming_csv_AI_points = pd.DataFrame({
        'transformed_fpr_normal': transformed_fpr_normal_AI,
        'transformed_tpr_normal': transformed_tpr_normal_AI
    })
    deming_csv_AI_points.to_csv(f'outputs/outputs_ID_{ID}_density_{density}/deming_AI_points.csv', index=False)

    ### Fit mu and sigma for non-diseased and diseased normal distributions for normal AI ROC curve

    def generate_roc_curve(fp_distribution, tp_distribution, size=100):
        min_threshold = min(fp_distribution.ppf(0.01), tp_distribution.ppf(0.01))
        max_threshold = max(fp_distribution.ppf(0.99), tp_distribution.ppf(0.99))
        min_threshold_fixed = -5
        max_threshold_fixed = 5
        thresholds = np.linspace(min_threshold, max_threshold, size)
        thresholds_fixed = np.linspace(min_threshold_fixed, max_threshold_fixed, size)
        fpr = np.array([1 - fp_distribution.cdf(t) for t in thresholds])
        tpr = np.array([1 - tp_distribution.cdf(t) for t in thresholds])
        fpr_fixed = np.array([1 - fp_distribution.cdf(t) for t in thresholds_fixed])
        tpr_fixed = np.array([1 - tp_distribution.cdf(t) for t in thresholds_fixed])
        return fpr, tpr, fpr_fixed, tpr_fixed

    # Transform Deming regression lines back to ROC values
    F_nd = norm(loc=0, scale=1)
    F_d_AI = norm(loc=a_normal_AI / b_normal_AI, scale=1 / b_normal_AI)
    fpr_AI_fit, tpr_AI_fit, fpr_AI_fit_fixed, tpr_AI_fit_fixed = generate_roc_curve(F_nd, F_d_AI)

    plt.figure(figsize=(6, 6))
    plt.plot(fpr_AI_fit, tpr_AI_fit, color='green', label='Fitted ROC')
    plt.plot(empirical_fpr_AI, empirical_tpr_AI, color='red', label='Empirical ROC', lw=2)
    plt.xlabel('FPR')
    plt.ylabel('TPR')
    plt.title('Fitted Smooth ROC for AI')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'outputs/outputs_ID_{ID}_density_{density}/empirical_fitted_AI.png')
    # plt.show()

    r_AI = a_normal_AI / (1 - b_normal_AI)
    print('mu:', a_normal_AI / b_normal_AI, 'sigma:', 1 / b_normal_AI, 'mean-to-sigma ratio:', a_normal_AI / (1 - b_normal_AI))

    AI_FPR_TPR_fit = pd.DataFrame({
        'fpr_AI_fit': fpr_AI_fit,
        'tpr_AI_fit': tpr_AI_fit,
        'fpr_AI_fit_fixed': fpr_AI_fit_fixed,
        'tpr_AI_fit_fixed': tpr_AI_fit_fixed
    })
    AI_FPR_TPR_fit.to_csv(f'outputs/outputs_ID_{ID}_density_{density}/AI_FPR_TPR_fit.csv', index=False)

    ### If 'radscores' exists (not just BIRADS), transform empirical radiologist ROC curve using Deming regression

    if 'radscores' in df.columns:
        empirical_fpr_rad, empirical_tpr_rad, thresholds_binary = roc_curve(df['label'], df['radscores'])

        max_first_index = max(np.argmax(empirical_fpr_rad > 0), np.argmax(empirical_tpr_rad > 0))
        min_last_index = min(np.argmax(empirical_fpr_rad >= 1), np.argmax(empirical_tpr_rad >= 1))
        if max_first_index >= min_last_index:
            min_last_index = max_first_index + 1

        transformed_fpr_normal_rad = transform_roc_values(empirical_fpr_rad[max_first_index:min_last_index], 'normal')
        transformed_tpr_normal_rad = transform_roc_values(empirical_tpr_rad[max_first_index:min_last_index], 'normal')

        a_normal_rad, b_normal_rad = fit_deming_regression(transformed_fpr_normal_rad, transformed_tpr_normal_rad)

        min_transformed_fpr_rad = np.min(transformed_fpr_normal_rad)
        max_transformed_fpr_rad = np.max(transformed_fpr_normal_rad)
        x_values_normal_rad = np.linspace(min_transformed_fpr_rad - 0.1, max_transformed_fpr_rad + 0.1, 400)
        y_values_normal_rad = deming_model([a_normal_rad, b_normal_rad], x_values_normal_rad)

        plt.figure(figsize=(6, 6))
        plt.plot(x_values_normal_rad, y_values_normal_rad, color='blue', label='Deming Regression Line (Normal)')
        plt.scatter(transformed_fpr_normal_rad, transformed_tpr_normal_rad, color='blue')
        plt.xlabel('Transformed FPR')
        plt.ylabel('Transformed TPR')
        plt.title('Deming Regression for Radiologist')
        plt.legend()
        plt.grid(True)
        plt.savefig(f'outputs/outputs_ID_{ID}_density_{density}/Deming_rad.png')
        # plt.show()

        deming_csv_rad_line = pd.DataFrame({
            'x_values_normal': x_values_normal_rad,
            'y_values_normal': y_values_normal_rad
        })
        deming_csv_rad_line.to_csv(f'outputs/outputs_ID_{ID}_density_{density}/deming_rad_line.csv', index=False)

        deming_csv_rad_points = pd.DataFrame({
            'transformed_fpr_normal': transformed_fpr_normal_rad,
            'transformed_tpr_normal': transformed_tpr_normal_rad
        })
        deming_csv_rad_points.to_csv(f'outputs/outputs_ID_{ID}_density_{density}/deming_rad_points.csv', index=False)

    ### Fit mu and sigma for non-diseased and diseased normal distributions for normal radiologist ROC curve (using all radiologists).

    if 'radscores' in df.columns:
        F_nd = norm(loc=0, scale=1)
        F_d_rad = norm(loc=a_normal_rad / b_normal_rad, scale=1 / b_normal_rad)
        fpr_rad_fit, tpr_rad_fit, fpr_rad_fit_fixed, tpr_rad_fit_fixed = generate_roc_curve(F_nd, F_d_rad)

        plt.figure(figsize=(6, 6))
        plt.plot(fpr_rad_fit, tpr_rad_fit, color='green', label='Fitted ROC')
        plt.plot(empirical_fpr_rad, empirical_tpr_rad, color='red', label='Empirical ROC', lw=2)
        plt.xlabel('FPR')
        plt.ylabel('TPR')
        plt.title('Fitted Smooth ROC for Radiologist')
        plt.legend()
        plt.grid(True)
        plt.savefig(f'outputs/outputs_ID_{ID}_density_{density}/empirical_fitted_rad.png')
        # plt.show()

        r_rad = a_normal_rad / (1 - b_normal_rad)

        rad_FPR_TPR_fit_continuous = pd.DataFrame({
            'fpr_rad_fit': fpr_rad_fit,
            'tpr_rad_fit': tpr_rad_fit,
            'fpr_rad_fit_fixed': fpr_rad_fit_fixed,
            'tpr_rad_fit_fixed': tpr_rad_fit_fixed
        })
        rad_FPR_TPR_fit_continuous.to_csv(f'outputs/outputs_ID_{ID}_density_{density}/rad_FPR_TPR_fit_continuous.csv', index=False)

    ### If only have BI-RADS scores and not radiologist scores on a finer scale, fit radiologist ROC curve to single point.
    ### Here we assume the same mean-to-sigma ratio "r" value for the radiologist curve as for the AI curve.

    if 'BIRADS' in df.columns:
        empirical_fpr_rad_birads, empirical_tpr_rad_birads = cal_birads_point(df)

        # Define the objective function: we want the ROC curve to pass through the radiologist (FPF, TPF) point.
        def objective_function(a, x, y, r):
            loc = a / (1 - a / r)
            scale = 1 / (1 - a / r)
            return abs(norm.sf(x, loc=loc, scale=scale) - y)

        result = minimize_scalar(objective_function, args=(norm.ppf(1 - empirical_fpr_rad_birads),
                                                             empirical_tpr_rad_birads, r_AI))
        if result.success:
            a_solution = result.x
            print("Solution for 'a':", a_solution)
        else:
            print("Failed to find a solution.")
            return

        # Generate radiologist ROC curve using the updated parameters and assuming same mean-to-sigma ratio as fitted AI curve (r_AI)
        F_nd_rad_birads = norm(loc=0, scale=1)
        F_d_rad_birads = norm(loc=a_solution / (1 - a_solution / r_AI), scale=1 / (1 - a_solution / r_AI))
        fpr_rad_fit_birads, tpr_rad_fit_birads, fpr_rad_fit_birads_fixed, tpr_rad_fit_birads_fixed = generate_roc_curve(F_nd_rad_birads, F_d_rad_birads)

        # Compute AUC under the fitted ROC curve
        AUC = np.abs(np.trapz(tpr_rad_fit_birads, fpr_rad_fit_birads))
        print("AUC under the fitted radiologist ROC curve:", AUC)

        # Plot ROC curves
        plt.figure(figsize=(6, 6))
        plt.plot(fpr_rad_fit_birads, tpr_rad_fit_birads, color='green', label='Smooth Normal Radiologist ROC curve')
        plt.plot(fpr_AI_fit, tpr_AI_fit, label='Smooth Normal AI ROC curve')
        plt.scatter(empirical_fpr_rad_birads, empirical_tpr_rad_birads, color='red', label='Radiologist (FPF,TPF)')
        plt.xlabel('FPF')
        plt.ylabel('TPF')
        plt.title('Normal ROC Curve with Fitted Radiologist FPF/TPF')
        plt.grid(True)
        plt.legend()
        plt.savefig(f'outputs/outputs_ID_{ID}_density_{density}/smooth_AI_rad.png')
        # plt.show()

        rad_FPR_TPR_fit_birads = pd.DataFrame({
            'rad_fpr_fit_birads': fpr_rad_fit_birads,
            'rad_tpr_fit_birads': tpr_rad_fit_birads,
            'rad_fpr_fit_birads_fixed': fpr_rad_fit_birads_fixed,
            'rad_tpr_fit_birads_fixed': tpr_rad_fit_birads_fixed
        })
        rad_FPR_TPR_fit_birads.to_csv(f'outputs/outputs_ID_{ID}_density_{density}/rad_FPR_TPR_fit_birads.csv', index=False)

    # Save a, b values for radiologist and AI ROC curves so that they can be reproduced.
    if 'radscores' in df.columns and 'BIRADS' not in df.columns:
        a_b_values = pd.DataFrame({
            "a_normal_AI": [a_normal_AI],
            "a_normal_rad": [a_normal_rad],
            "b_normal_AI": [b_normal_AI],
            "b_normal_rad": [b_normal_rad]
        })
    elif 'radscores' in df.columns and 'BIRADS' in df.columns:
        a_b_values = pd.DataFrame({
            "a_normal_AI": [a_normal_AI],
            "a_normal_rad": [a_normal_rad],
            "b_normal_AI": [b_normal_AI],
            "b_normal_rad": [b_normal_rad],
            "a_solution_BIRADS": [a_solution]
        })
    elif 'radscores' not in df.columns and 'BIRADS' in df.columns:
        a_b_values = pd.DataFrame({
            "a_normal_AI": [a_normal_AI],
            "b_normal_AI": [b_normal_AI],
            "a_solution_BIRADS": [a_solution]
        })

    a_b_values.to_csv(f'outputs/outputs_ID_{ID}_density_{density}/a_b_values.csv', index=False)

    ### Compute Correlations between BIRADS and AI scores

    def correlations_BIRADS(dfgroup):
        if 'radID' in dfgroup.columns and condition_reader == 'False':
            if are_all_same(dfgroup, 'radID'):
                groupName = dfgroup['radID'].iloc[0]
            else:
                groupName = 'All'
        else:
            groupName = 'All'

        if 'BIRADS' in dfgroup.columns:
            nd_rows = dfgroup[dfgroup['label'] == 0]
            d_rows = dfgroup[dfgroup['label'] == 1]

            F_nd = norm(loc=0, scale=1)
            F_d_AI = norm(loc=a_normal_AI / b_normal_AI, scale=1 / b_normal_AI)
            F_d_rad_birads = norm(loc=a_solution / (1 - a_solution / r_AI), scale=1 / (1 - a_solution / r_AI))

            rad_labels_nd_binary = nd_rows['rad']
            rad_labels_nd = nd_rows['rad_ANB']
            continuous_scores_nd = nd_rows['AIscores']
            rad_labels_d_binary = d_rows['rad']
            rad_labels_d = d_rows['rad_ANB']
            continuous_scores_d = d_rows['AIscores']
            all_labels_binary = dfgroup['rad']
            all_labels = dfgroup['rad_ANB']
            all_scores = dfgroup['AIscores']

            q_rad_labels_nd_binary = norm.ppf(F_nd.cdf(nd_rows['rad']))
            q_rad_labels_nd = norm.ppf(F_nd.cdf(nd_rows['rad_ANB']))
            q_continuous_scores_nd = norm.ppf(F_nd.cdf(nd_rows['AIscores']))
            q_rad_labels_d_binary = norm.ppf(F_d_rad_birads.cdf(d_rows['rad']))
            q_rad_labels_d = norm.ppf(F_d_rad_birads.cdf(d_rows['rad_ANB']))
            q_continuous_scores_d = norm.ppf(F_d_AI.cdf(d_rows['AIscores']))
            q_rad_labels_binary_all = np.hstack((q_rad_labels_nd_binary, q_rad_labels_d_binary))
            q_rad_labels_all = np.hstack((q_rad_labels_nd, q_rad_labels_d))
            q_continuous_scores_all = np.hstack((continuous_scores_nd, q_continuous_scores_d))

            mask_nd = (~np.isnan(q_continuous_scores_nd)) & (~np.isnan(q_rad_labels_nd_binary)) & \
                      np.isfinite(q_continuous_scores_nd) & np.isfinite(q_rad_labels_nd_binary)
            mask_d = (~np.isnan(q_continuous_scores_d)) & (~np.isnan(q_rad_labels_d_binary)) & \
                     np.isfinite(q_continuous_scores_d) & np.isfinite(q_rad_labels_d_binary)
            mask_all = (~np.isnan(q_continuous_scores_all)) & (~np.isnan(q_rad_labels_all)) & \
                       np.isfinite(q_continuous_scores_all) & np.isfinite(q_rad_labels_all)

            q_rho_all_binary, p_value_q_rho_all_binary = pearsonr(q_continuous_scores_all[mask_all], q_rad_labels_binary_all[mask_all])

            if len(nd_rows) > 1:
                rho_nd_binary, p_value_rho_nd_binary = pearsonr(continuous_scores_nd, rad_labels_nd_binary)
                rho_nd, p_value_rho_nd = pearsonr(continuous_scores_nd, rad_labels_nd)
                tau_nd_binary, p_value_nd_binary = kendalltau(continuous_scores_nd, rad_labels_nd_binary)
                tau_nd, p_value_nd = kendalltau(continuous_scores_nd, rad_labels_nd)
                q_rho_nd_binary, p_value_q_rho_nd_binary = pearsonr(q_continuous_scores_nd[mask_nd], q_rad_labels_nd_binary[mask_nd])
                q_rho_nd, p_value_q_rho_nd = pearsonr(q_continuous_scores_nd[mask_nd], q_rad_labels_nd[mask_nd])
            else:
                rho_nd_binary, p_value_rho_nd_binary = 1, 0
                rho_nd, p_value_rho_nd = 1, 0
                tau_nd_binary, p_value_nd_binary = 1, 0
                tau_nd, p_value_nd = 1, 0
                q_rho_nd_binary, p_value_q_rho_nd_binary = 1, 0
                q_rho_nd, p_value_q_rho_nd = 1, 0

            if len(d_rows) > 1:
                rho_d_binary, p_value_rho_d_binary = pearsonr(continuous_scores_d, rad_labels_d_binary)
                rho_d, p_value_rho_d = pearsonr(continuous_scores_d, rad_labels_d)
                tau_d_binary, p_value_d_binary = kendalltau(continuous_scores_d, rad_labels_d_binary)
                tau_d, p_value_d = kendalltau(continuous_scores_d, rad_labels_d)
                q_rho_d_binary, p_value_q_rho_d_binary = pearsonr(q_continuous_scores_d[mask_d], q_rad_labels_d_binary[mask_d])
                q_rho_d, p_value_q_rho_d = pearsonr(q_continuous_scores_d[mask_d], q_rad_labels_d[mask_d])
            else:
                rho_d_binary, p_value_rho_d_binary = 1, 0
                rho_d, p_value_rho_d = 1, 0
                tau_d_binary, p_value_d_binary = 1, 0
                tau_d, p_value_d = 1, 0
                q_rho_d_binary, p_value_q_rho_d_binary = 1, 0
                q_rho_d, p_value_q_rho_d = 1, 0

            if len(nd_rows) + len(d_rows) > 1:
                rho_all_binary, p_value_rho_all_binary = pearsonr(all_scores, all_labels_binary)
                rho_all, p_value_rho_all = pearsonr(all_scores, all_labels)
                tau_all_binary, p_value_all_binary = kendalltau(all_scores, all_labels_binary)
                tau_all, p_value_all = kendalltau(all_scores, all_labels)
                q_rho_all_binary, p_value_q_rho_all_binary = pearsonr(q_continuous_scores_all[mask_all], q_rad_labels_binary_all[mask_all])
                q_rho_all, p_value_q_rho_all = pearsonr(q_continuous_scores_all[mask_all], q_rad_labels_all[mask_all])
            else:
                rho_all_binary, p_value_rho_all_binary = 1, 0
                rho_all, p_value_rho_all = 1, 0
                tau_all_binary, p_value_all_binary = 1, 0
                tau_all, p_value_all = 1, 0
                q_rho_all_binary, p_value_q_rho_all_binary = 1, 0
                q_rho_all, p_value_q_rho_all = 1, 0

            data = {
                'Statistic': ["Pearson's rho"] * 6 + ["Kendall's tau"] * 6 + ["Pearson's rho, inv normal"] * 6,
                'Group': ['nd', 'd', 'all'] * 6,
                'Value': [rho_nd_binary, rho_d_binary, rho_all_binary, rho_nd, rho_d, rho_all,
                          tau_nd_binary, tau_d_binary, tau_all_binary, tau_nd, tau_d, tau_all,
                          q_rho_nd_binary, q_rho_d_binary, q_rho_all_binary, q_rho_nd, q_rho_d, q_rho_all],
                'P-value': [p_value_rho_nd_binary, p_value_rho_d_binary, p_value_rho_all_binary, p_value_rho_nd,
                            p_value_rho_d, p_value_rho_all, p_value_nd_binary, p_value_d_binary, p_value_all_binary,
                            p_value_nd, p_value_d, p_value_all, p_value_q_rho_nd_binary, p_value_q_rho_d_binary,
                            p_value_q_rho_all_binary, p_value_q_rho_nd, p_value_q_rho_d, p_value_q_rho_all]
            }

            correlation_results_BIRADS = pd.DataFrame(data)
            correlation_results_BIRADS_str = correlation_results_BIRADS.to_markdown(index=False)
            print(f'READER GROUP: {groupName}')
            print(correlation_results_BIRADS_str)

            correlation_results_BIRADS.to_csv(f'outputs/outputs_ID_{ID}_density_{density}/correlation_results_BIRADS_{groupName}.csv', index=False)

    correlations_BIRADS(df)

    if 'radID' in df.columns and 'BIRADS' in df.columns:
        random_IDs = df.groupby('radID').ngroup().add(1)  # Replace radiologist IDs with random numbers
        df['radID'] = random_IDs
        grouped = df.groupby('radID')
        grouped.apply(correlations_BIRADS)

    ### Compute Correlations between finer-scale radiologist scores and AI scores

    def correlations_continuous(dfgroup):
        if 'radID' in dfgroup.columns and condition_reader == 'False':
            if are_all_same(dfgroup, 'radID'):
                groupName = dfgroup['radID'].iloc[0]
            else:
                groupName = 'All'
        else:
            groupName = 'All'

        if 'radscores' in df.columns:
            F_nd = norm(loc=0, scale=1)
            F_d_AI = norm(loc=a_normal_AI / b_normal_AI, scale=1 / b_normal_AI)
            F_d_rad = norm(loc=a_normal_rad / b_normal_rad, scale=1 / b_normal_rad)

            nd_rows_rad = df[df['label'] == 0]
            d_rows_rad = df[df['label'] == 1]

            rad_scores_nd = nd_rows_rad['radscores']
            rad_scores_d = d_rows_rad['radscores']
            AI_scores_nd = nd_rows_rad['AIscores']
            AI_scores_d = d_rows_rad['AIscores']
            all_scores_rad = df['radscores']
            all_scores_AI = df['AIscores']

            q_rad_labels_nd = norm.ppf(F_nd.cdf(nd_rows_rad['radscores']))
            q_continuous_scores_nd = norm.ppf(F_nd.cdf(nd_rows_rad['AIscores']))
            q_rad_labels_d = norm.ppf(F_d_rad.cdf(d_rows_rad['radscores']))
            q_continuous_scores_d = norm.ppf(F_d_AI.cdf(d_rows_rad['AIscores']))
            q_rad_labels_all = np.hstack((q_rad_labels_nd, q_rad_labels_d))
            q_continuous_scores_all = np.hstack((AI_scores_nd, q_continuous_scores_d))

            if len(nd_rows_rad) > 1:
                rho_nd_continuous, p_value_rho_nd_continuous = pearsonr(rad_scores_nd, AI_scores_nd)
                tau_nd_continuous, p_value_tau_nd_continuous = kendalltau(rad_scores_nd, AI_scores_nd)
                q_rho_nd, p_value_q_rho_nd = pearsonr(q_continuous_scores_nd, q_rad_labels_nd)
            else:
                rho_nd_continuous, p_value_rho_nd_continuous = 1, 0
                tau_nd_continuous, p_value_tau_nd_continuous = 1, 0
                q_rho_nd, p_value_q_rho_nd = 1, 0

            if len(d_rows_rad) > 1:
                rho_d_continuous, p_value_rho_d_continuous = pearsonr(rad_scores_d, AI_scores_d)
                tau_d_continuous, p_value_tau_d_continuous = kendalltau(rad_scores_d, AI_scores_d)
                q_rho_d, p_value_q_rho_d = pearsonr(q_continuous_scores_d, q_rad_labels_d)
            else:
                rho_d_continuous, p_value_rho_d_continuous = 1, 0
                tau_d_continuous, p_value_tau_d_continuous = 1, 0
                q_rho_d, p_value_q_rho_d = 1, 0

            if len(nd_rows_rad) + len(d_rows_rad) > 1:
                rho_all_continuous, p_value_rho_all_continuous = pearsonr(all_scores_rad, all_scores_AI)
                tau_all_continuous, p_value_tau_all_continuous = kendalltau(all_scores_rad, all_scores_AI)
                q_rho_all, p_value_q_rho_all = pearsonr(q_continuous_scores_all, q_rad_labels_all)
            else:
                rho_all_continuous, p_value_rho_all_continuous = 1, 0
                tau_all_continuous, p_value_tau_all_continuous = 1, 0
                q_rho_all, p_value_q_rho_all = 1, 0

            data = {
                'Statistic': ["Pearson's rho"] * 3 + ["Kendall's tau"] * 3 + ["Pearson's rho, inv normal"] * 3,
                'Group': ['nd', 'd', 'all'] * 3,
                'Value': [rho_nd_continuous, rho_d_continuous, rho_all_continuous,
                          tau_nd_continuous, tau_d_continuous, tau_all_continuous,
                          q_rho_nd, q_rho_d, q_rho_all],
                'P-value': [p_value_rho_nd_continuous, p_value_rho_d_continuous, p_value_rho_all_continuous,
                            p_value_tau_nd_continuous, p_value_tau_d_continuous, p_value_tau_all_continuous,
                            p_value_q_rho_nd, p_value_q_rho_d, p_value_q_rho_all]
            }

            correlation_results_continuous = pd.DataFrame(data)
            correlation_results_continuous.to_csv(f'outputs/outputs_ID_{ID}_density_{density}/correlation_results_continuous_{groupName}.csv', index=False)

            print(f'READER GROUP: {groupName}')
            correlation_results_continuous_str = correlation_results_continuous.to_markdown(index=False)
            print(correlation_results_continuous_str)

    correlations_continuous(df)

    if 'radID' in df.columns and 'radscores' in df.columns:
        random_IDs = df.groupby('radID').ngroup().add(1)  # Replace radiologist IDs with random numbers
        df['radID'] = random_IDs
        grouped = df.groupby('radID')
        grouped.apply(correlations_continuous)

    ### Compute correlations between radiologists (if applicable)

    if 'imageID' in df.columns and 'radID' in df.columns:
        def reader_scores(df, type):
            if len(df['label']) == len(df[df['label'] == 1]):
                name = 'diseased'
            elif len(df['label']) == len(df[df['label'] == 0]):
                name = 'non-diseased'
            else:
                name = 'all'
            pairs = list(combinations(df['radID'].unique(), 2))
            correlations = {}
            taucorrelations = {}
            for pair in pairs:
                rad1, rad2 = pair
                scores_rad1 = []
                scores_rad2 = []
                for image in df['imageID'].unique():
                    subset = df[df['imageID'] == image]
                    if type == 'birads':
                        score_rad1 = subset.loc[subset['radID'] == rad1, 'rad'].values
                        score_rad2 = subset.loc[subset['radID'] == rad2, 'rad'].values
                    if type == 'continuous':
                        score_rad1 = subset.loc[subset['radID'] == rad1, 'radscores'].values
                        score_rad2 = subset.loc[subset['radID'] == rad2, 'radscores'].values
                    if len(score_rad1) > 0 and len(score_rad2) > 0:
                        scores_rad1.append(score_rad1[0])
                        scores_rad2.append(score_rad2[0])
                if len(scores_rad1) > 1 and len(scores_rad2) > 1:
                    correlation, _ = pearsonr(scores_rad1, scores_rad2)
                    taucorrelation, _ = kendalltau(scores_rad1, scores_rad2)
                    correlations[pair] = correlation
                    taucorrelations[pair] = taucorrelation
            for pair, correlation in correlations.items():
                print(f"Pearson {name} correlation between radiologists {pair[0]} and {pair[1]}: {correlation}, {type}")
            for pair, taucorrelation in taucorrelations.items():
                print(f"Kendall's tau {name} correlation between radiologists {pair[0]} and {pair[1]}: {taucorrelation}, {type}")
            results = pd.DataFrame({
                'Pair': [f"{pair[0]}-{pair[1]}" for pair in correlations.keys()],
                'Pearson Correlation': list(correlations.values()),
                'Kendall Tau Correlation': list(taucorrelations.values())
            })
            return results

        if 'BIRADS' in df.columns:
            all_results_birads = reader_scores(df, 'birads')
            print(' ')
            diseased_results_birads = reader_scores(df[df['label'] == 1], 'birads')
            print(' ')
            non_diseased_results_birads = reader_scores(df[df['label'] == 0], 'birads')
            print(' ')

            all_results_birads.to_csv(f'outputs/outputs_ID_{ID}_density_{density}/all_correlation_results_birads.csv', index=False)
            diseased_results_birads.to_csv(f'outputs/outputs_ID_{ID}_density_{density}/diseased_correlation_results_birads.csv', index=False)
            non_diseased_results_birads.to_csv(f'outputs/outputs_ID_{ID}_density_{density}/non_diseased_correlation_results_birads.csv', index=False)

        if 'radscores' in df.columns:
            all_results_continuous = reader_scores(df, 'continuous')
            print(' ')
            diseased_results_continuous = reader_scores(df[df['label'] == 1], 'continuous')
            print(' ')
            non_diseased_results_continuous = reader_scores(df[df['label'] == 0], 'continuous')

            all_results_continuous.to_csv(f'outputs/outputs_ID_{ID}_density_{density}/all_correlation_results_continuous.csv', index=False)
            diseased_results_continuous.to_csv(f'outputs/outputs_ID_{ID}_density_{density}/diseased_correlation_results_continuous.csv', index=False)
            non_diseased_results_continuous.to_csv(f'outputs/outputs_ID_{ID}_density_{density}/non_diseased_correlation_results_continuous.csv', index=False)

    print('Warning: the remaining code generates a csv file with radiologist performance under many rule-out/rule-in scenarios and outputs a few illustrative plots. This may take several minutes to run depending on the size of your dataset.')

    ### Generate a csv file with combinations of rule-out and rule-in AI FPF cutoffs and corresponding empirical radiologist performance under that rule-out+rule-in scenario.
    ### Plot a few of the results along with theoretical rule-out+rule-in curves from copulas.

    if 'BIRADS' in df.columns:
        min_AI_score = df['AIscores'].min()
        max_AI_score = df['AIscores'].max()
        rule_out_cutoffs = beta.ppf(np.linspace(min_AI_score, max_AI_score, 100), 2, 2)
        rule_in_cutoffs = beta.ppf(np.linspace(min_AI_score, max_AI_score, 100), 2, 2)
        rule_out_cutoff = []
        rule_in_cutoff = []
        joint_rad_FPF_list = []
        joint_rad_TPF_list = []
        AI_FPF_ruleout_list = []
        AI_FPF_rulein_list = []
        perc_ruledout_list = []
        perc_ruledin_list = []

        for x in rule_out_cutoffs:
            for y in rule_in_cutoffs:
                if x < y:
                    df['jointrad'] = 0
                    AI_FPF_ruleout = df[(df['label'] == 0) & (df['AIscores'] > x)].shape[0] / (
                        df[(df['label'] == 0) & (df['AIscores'] > x)].shape[0] + df[(df['label'] == 0) & (df['AIscores'] < x)].shape[0])
                    AI_FPF_rulein = df[(df['label'] == 0) & (df['AIscores'] > y)].shape[0] / (
                        df[(df['label'] == 0) & (df['AIscores'] > y)].shape[0] + df[(df['label'] == 0) & (df['AIscores'] < y)].shape[0])
                    perc_ruledout = df[(df['AIscores'] < x)].shape[0] / len(df)
                    perc_ruledin = df[(df['AIscores'] > y)].shape[0] / len(df)
                    df.loc[(df['AIscores'] > x) & (df['rad'] == 1), 'jointrad'] = 1
                    df.loc[(df['AIscores'] > y), 'jointrad'] = 1
                    TP = df[(df['jointrad'] == 1) & (df['label'] == 1)].shape[0]
                    FP = df[(df['jointrad'] == 1) & (df['label'] == 0)].shape[0]
                    TN = df[(df['jointrad'] == 0) & (df['label'] == 0)].shape[0]
                    FN = df[(df['jointrad'] == 0) & (df['label'] == 1)].shape[0]
                    joint_rad_FPF = FP / (FP + TN)
                    joint_rad_TPF = TP / (TP + FN)

                    joint_rad_FPF_list.append(joint_rad_FPF)
                    joint_rad_TPF_list.append(joint_rad_TPF)
                    AI_FPF_ruleout_list.append(AI_FPF_ruleout)
                    AI_FPF_rulein_list.append(AI_FPF_rulein)
                    perc_ruledout_list.append(perc_ruledout)
                    perc_ruledin_list.append(perc_ruledin)
                    rule_out_cutoff.append(x)
                    rule_in_cutoff.append(y)

        result_df_birads = pd.DataFrame({
            'joint_rad_FPF': joint_rad_FPF_list,
            'joint_rad_TPF': joint_rad_TPF_list,
            'AI_FPF_ruleout': AI_FPF_ruleout_list,
            'AI_FPF_rulein': AI_FPF_rulein_list,
            'perc_ruledout': perc_ruledout_list,
            'perc_ruledin': perc_ruledin_list,
            'rule_out_list': rule_out_cutoff,
            'rule_in_list': rule_in_cutoff
        })

        result_df_birads.to_csv(f'outputs/outputs_ID_{ID}_density_{density}/AI_FPF_rad_performance_birads.csv', index=False)

    ### Generate a csv file with empirical joint ROC curves for the radiologist under different rule-out+rule-in scenarios.

    if 'radscores' in df.columns:
        min_AI_score = df['AIscores'].min()
        max_AI_score = df['AIscores'].max()
        min_rad_score = df['radscores'].min()
        max_rad_score = df['radscores'].max()
        lowest_score = df['radscores'].min() - 1
        highest_score = df['radscores'].max() + 1

        rule_out_cutoffs = beta.ppf(np.linspace(min_AI_score, max_AI_score, 100), 2, 2)
        rule_in_cutoffs = beta.ppf(np.linspace(min_AI_score, max_AI_score, 100), 2, 2)
        rule_out_cutoff = []
        rule_in_cutoff = []
        joint_rad_FPF_list = []
        joint_rad_TPF_list = []
        AI_FPF_ruleout_list = []
        AI_FPF_rulein_list = []
        perc_ruledout_list = []
        perc_ruledin_list = []

        for x in rule_out_cutoffs:
            for y in rule_in_cutoffs:
                if x < y:
                    dfnew = df.copy()
                    dfnew['jointrad'] = df['radscores']
                    AI_FPF_ruleout = df[(df['label'] == 0) & (df['AIscores'] > x)].shape[0] / (
                        df[(df['label'] == 0) & (df['AIscores'] > x)].shape[0] + df[(df['label'] == 0) & (df['AIscores'] < x)].shape[0])
                    AI_FPF_rulein = df[(df['label'] == 0) & (df['AIscores'] > y)].shape[0] / (
                        df[(df['label'] == 0) & (df['AIscores'] > y)].shape[0] + df[(df['label'] == 0) & (df['AIscores'] < y)].shape[0])
                    perc_ruledout = df[(df['AIscores'] < x)].shape[0] / len(df)
                    perc_ruledin = df[(df['AIscores'] > y)].shape[0] / len(df)
                    dfnew.loc[df['AIscores'] < x, 'jointrad'] = lowest_score
                    dfnew.loc[df['AIscores'] > y, 'jointrad'] = highest_score
                    empirical_fpr_rad_joint, empirical_tpr_rad_joint, thresholds_binary = roc_curve(dfnew['label'], dfnew['jointrad'])
                    empirical_fpr_rad_joint_list = list(empirical_fpr_rad_joint)
                    empirical_tpr_rad_joint_list = list(empirical_tpr_rad_joint)

                    AI_FPF_ruleout_list.append(AI_FPF_ruleout)
                    AI_FPF_rulein_list.append(AI_FPF_rulein)
                    perc_ruledout_list.append(perc_ruledout)
                    perc_ruledin_list.append(perc_ruledin)
                    rule_out_cutoff.append(x)
                    rule_in_cutoff.append(y)
                    joint_rad_FPF_list.append(empirical_fpr_rad_joint_list)
                    joint_rad_TPF_list.append(empirical_tpr_rad_joint_list)

        result_df_continuous = pd.DataFrame({
            'joint_rad_FPF': joint_rad_FPF_list,
            'joint_rad_TPF': joint_rad_TPF_list,
            'AI_FPF_ruleout': AI_FPF_ruleout_list,
            'AI_FPF_rulein': AI_FPF_rulein_list,
            'perc_ruledout': perc_ruledout_list,
            'perc_ruledin': perc_ruledin_list,
            'rule_out_list': rule_out_cutoff,
            'rule_in_list': rule_in_cutoff
        })

        result_df_continuous.to_csv(f'outputs/outputs_ID_{ID}_density_{density}/AI_FPF_rad_performance_continuous.csv', index=False)


if 'radID' in df_orig.columns:
    rad_set = ['True', 'False']
else:
    rad_set = ['False']

if 'density' in df_orig.columns:
    dens_set = ['True', 'False']
else:
    dens_set = ['False']

for condition_reader in rad_set:
    for condition_density in dens_set:
        if condition_reader == 'True' and 'radID' in df_orig.columns and condition_density == 'False':
            random_IDs = df_orig.groupby('radID').ngroup().add(1)  # Replace radiologist IDs with random numbers
            df_orig['radID'] = random_IDs
            for ID in random_IDs.unique():
                print('ID:', ID)
                if not os.path.exists(f'outputs/outputs_ID_{ID}_density_full/'):
                    os.mkdir(f'outputs/outputs_ID_{ID}_density_full/')
                df = df_orig[df_orig['radID'] == ID]
                full_analysis(df, ID, 'full')

        elif condition_reader == 'False' and condition_density == 'True' and 'density' in df_orig.columns:
            for dens in df_orig['density'].unique():
                print('Density:', dens)
                if not os.path.exists(f'outputs/outputs_ID_full_density_{dens}/'):
                    os.mkdir(f'outputs/outputs_ID_full_density_{dens}/')
                df = df_orig[df_orig['density'] == dens]
                full_analysis(df, 'full', dens)

        elif condition_reader == 'True' and 'radID' in df_orig.columns and condition_density == 'True' and 'density' in df_orig.columns:
            random_IDs = df_orig.groupby('radID').ngroup().add(1)  # Replace radiologist IDs with random numbers
            df_orig['radID'] = random_IDs
            for ID in random_IDs.unique():
                for dens in df_orig['density'].unique():
                    print('ID:', ID)
                    print('Density:', dens)
                    if not os.path.exists(f'outputs/outputs_ID_{ID}_density_{dens}/'):
                        os.mkdir(f'outputs/outputs_ID_{ID}_density_{dens}/')
                    df_new = df_orig[df_orig['radID'] == ID]
                    df = df_new[df_new['density'] == dens]
                    full_analysis(df, ID, dens)

        else:
            if not os.path.exists(f'outputs/outputs_ID_full_density_full/'):
                os.mkdir(f'outputs/outputs_ID_full_density_full/')
            full_analysis(df_orig, 'full', 'full')
