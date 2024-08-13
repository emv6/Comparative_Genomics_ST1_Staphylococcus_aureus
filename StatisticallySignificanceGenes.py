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

#Iterating over the dataframe 
    for index, row in df.iterrows():
        print(f"Processing Gene: {row['Gene']}")

        #Contingency Table 
        table = np.array([
            [row['Bovine_Presence'], row['Bovine_Absence']],
            [row['Human_Presence'], row['Human_Absence']]
        ])

        try:
            #Chi-Square 
            chi2_stat, p_value_chi2, _, expected = stats.chi2_contingency(table, correction=False)
            #If any frequency is less than 5 
            if np.any(expected < 5):
                test_name = 'Fisher\'s exact'
                oddsratio_fisher, p_value_fisher = stats.fisher_exact(table, alternative ='two-sided')
                p_value = p_value_fisher
                chi2_stat = np.nan
            else:
                test_name = 'Chi-Square'
                p_value = p_value_chi2
                p_value_fisher = np.nan 

        #If there is any error result to Fisher's Exact Test 
        except Exception as e: 
            print(f"Exception encountered: {e}")
            test_name = 'Fisher\'s Exact'
            oddsratio_fisher, p_value_fisher = stats.fisher_exact(table, alternative ='two-sided')
            p_value = p_value_fisher 
            chi2_stat = np.nan 
            p_value_chi2 = np.nan 

        #Calculate Odds Ratio and 95% CI 
        a, b, c, d = table.flatten()
        odds_ratio, (ci_lower, ci_upper) = odds_ratio_and_ci(a, b, c, d)

         #If genes are significant say True otherwise False 
        is_significant = p_value < 0.05

    #Append the results for each gene
        results.append({
            'Gene': row['Gene'], 
            'Test Used': test_name, 
            'P-Value': p_value, 
            'Chi-Square Statistic': chi2_stat,
            'Chi-Square P value': p_value_chi2,
            'Fisher\'s Exact P-Value': p_value_fisher,
            'Odds Ratio': odds_ratio, 
            '95% CI Lower': ci_lower, 
            '95% CI Upper': ci_upper,
            'Significance': is_significant,
            'Expected': expected
        })

#Convert results to DataFrame 
    results_df = pd.DataFrame(results)

    return results_df

#Analyse data and print results 
results_df = analyze_genes(df)
print("Results DataFrame:")
print(results_df)
results_df.to_csv('AMR_Results.csv', index = False)