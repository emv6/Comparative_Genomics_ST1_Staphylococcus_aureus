#!/usr/bin/env python 3 

def print_header ():
    header = """
======================================
Script Name is Statistical_Significance_GenesMGE.py 
Description:
    This script performs statistical analysis on Gene Presence and Absence data 
    including Chi-Square, Fisher's Exact Test, Odds Ratio and 95% Confidence Interval for each gene. 
    Which statistical test to use is based on the expected frequencies. 

Author: Emma Voss
Email: emma.voss@lic.co.nz
Version 1.1

=======================================

"""
    print(header)
print_header()

import pandas as pd 
import numpy as np 
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

# Calculate Odds Ratio and 95% CI
def odds_ratio_and_ci(a, b, c, d, confidence_level=0.95):
    
    # Continuity Correction 0.5
    if a == 0 or b == 0 or c == 0 or d == 0:
        a += 0.5
        b += 0.5
        c += 0.5
        d += 0.5
    #Odds Ratio
    odds_ratio = (a * d) / (b * c) if (b * c) != 0 else np.nan
    
    # Standard error of the odds ratio
    se = np.sqrt((1/a) + (1/b) + (1/c) + (1/d))

    # Z-score for CI
    z = stats.norm.ppf(1 - (1 - confidence_level) / 2)

    # CI for odds ratio
    ci_lower = np.exp(np.log(odds_ratio) - z * se) if odds_ratio > 0 else np.nan
    ci_upper = np.exp(np.log(odds_ratio) + z * se) if odds_ratio > 0 else np.nan
    
    return odds_ratio, (ci_lower, ci_upper)

file_path = 'Contingency_Table.xlsx' 
sheet_name = 'Sheet1'
df = pd.read_excel(file_path, sheet_name=sheet_name)

print(df)

def analyze_genes(df):
    results = []
    raw_p_values = []

    for index, row in df.iterrows():
        gene = row['Gene']
        bovine_presence = row['Bovine_Presence']
        human_presence = row['Human_Presence']
        
        total_presence = bovine_presence + human_presence

        # Skip if total presence < 5
        if total_presence < 5:
            print(f"Skipping {gene} due to low total presence ({total_presence})")
            continue

        table = np.array([
            [bovine_presence, row['Bovine_Absence']],
            [human_presence, row['Human_Absence']]
        ])

        try:
            chi2_stat, p_value_chi2, _, expected = stats.chi2_contingency(table, correction=False)
            if np.any(expected < 5):
                test_name = 'Fisher\'s Exact'
                oddsratio_fisher, p_value_fisher = stats.fisher_exact(table)
                p_value = p_value_fisher
                chi2_stat = np.nan
            else:
                test_name = 'Chi-Square'
                p_value = p_value_chi2
                p_value_fisher = np.nan
        except Exception as e:
            print(f"Exception for {gene}: {e}")
            test_name = 'Fisher\'s Exact'
            oddsratio_fisher, p_value_fisher = stats.fisher_exact(table)
            p_value = p_value_fisher
            chi2_stat = np.nan
            p_value_chi2 = np.nan

        # Odds ratio and CI
        a, b, c, d = table.flatten()
        odds_ratio, (ci_lower, ci_upper) = odds_ratio_and_ci(a, b, c, d)

        raw_p_values.append(p_value)

        results.append({
            'Gene': gene,
            'Test Used': test_name,
            'P-Value': p_value,
            'Chi-Square Statistic': chi2_stat,
            'Chi-Square P-Value': p_value_chi2,
            'Fisher\'s Exact P-Value': p_value_fisher,
            'Odds Ratio': odds_ratio,
            '95% CI Lower': ci_lower,
            '95% CI Upper': ci_upper
        })

    # Convert to DataFrame
    results_df = pd.DataFrame(results)

    # Apply FDR correction
    if not results_df.empty:
        fdr_results = multipletests(raw_p_values, alpha=0.05, method='fdr_bh')
        results_df['FDR-adjusted P-Value'] = fdr_results[1]
        results_df['FDR Significant (<0.05)'] = fdr_results[0]
    else:
        results_df['FDR-adjusted P-Value'] = []
        results_df['FDR Significant (<0.05)'] = []

    return results_df

# Load Excel
file_path = 'Contingency_Table.xlsx'
sheet_name = 'Sheet1'
df = pd.read_excel(file_path, sheet_name=sheet_name)

# Analyze
results_df = analyze_genes(df)
print("Results DataFrame:")
print(results_df)

results_df.to_csv('AMR_Results.csv', index = False)
