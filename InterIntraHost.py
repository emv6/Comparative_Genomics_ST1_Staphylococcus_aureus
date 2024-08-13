import pandas as pd
import numpy as np
import scipy.stats as stats 
from itertools import combinations
import matplotlib.pyplot as plt
from statsmodels.stats.diagnostic import lilliefors

# Read the SNP distance matrix CSV file
snp_matrix_file = 'snpdist.csv' 
df = pd.read_csv(snp_matrix_file, index_col=0)

# Read the isolate-to-host mapping CSV file
isolate_host_mapping_file = 'Snp_Host.csv'
host_df = pd.read_csv(isolate_host_mapping_file)

# Create a dictionary mapping isolate labels to hosts
isolate_to_host = dict(zip(host_df['Isolate'], host_df['Host']))

# Debug print: Check the isolate_to_host dictionary
print("Isolate to Host Mapping Dictionary:")
print(isolate_to_host)

# Extract the labels
isolate_labels = df.columns.tolist()

# Convert the DataFrame to a NumPy array
snp_matrix = df.to_numpy()

# Function to compute average SNP distances
def compute_snp_distances(snp_matrix, isolate_labels, isolate_to_host):
   intra_distances = {}
   inter_distances = []

   # Get the set of unique hosts
   hosts = set(isolate_to_host.values())

   # Debug print: Check the set of hosts
   #print("Unique Hosts:")
   #print(hosts)

   # Initialize dictionaries for intra and inter distances: 
   for host in hosts:
       intra_distances[host] = []

   # Compute distances
   for (i, j) in combinations(range(len(snp_matrix)), 2):
       iso1, iso2 = isolate_labels[i], isolate_labels[j]
       host1, host2 = isolate_to_host.get(iso1), isolate_to_host.get(iso2)

       # Debug print: Check host values
       #print(f"Processing: {iso1} (Host: {host1}), {iso2} (Host: {host2})")

       # Skip if any host value is None (isolates not in the appropriate file)
       if host1 is None or host2 is None:
           continue
       distance = snp_matrix[i, j]
       if host1 == host2:
           intra_distances[host1].append(distance)
       else:
           inter_distances.append(distance)

   # Calculate averages
   avg_intra_distances = {host: np.mean(distances) if distances else float('nan')
                          for host, distances in intra_distances.items()}
   avg_inter_distance = np.mean(inter_distances) if inter_distances else float('nan')
   return intra_distances, inter_distances, avg_intra_distances, avg_inter_distance

# Calculate the SNP distances
intra_distances, inter_distances, avg_intra_distances, avg_inter_distance = compute_snp_distances(snp_matrix, isolate_labels, isolate_to_host)

# Check normality but can't use Shapiro-Wilk as first trial N was greater than 5000 so using Lilliefors
def check_normality(data, label):
   lilliefors_test = lilliefors(data)
   print(f"{label} - Lilliefors Test: D={lilliefors_test[0]}, p-value={lilliefors_test[1]}")
   plt.figure()
   stats.probplot(data, dist="norm", plot=plt)
   plt.title(f"Q-Q Plot for {label}")
   plt.show()
   return  lilliefors_test[1] > 0.05  # Return True if data is normally distributed
print("Normality Check:")
normality_results = {}

for host, distances in intra_distances.items(): 
    normality_results[host]= check_normality(distances, f"{host} Intra-Host Distances")
normality_results['Inter-host'] = check_normality(inter_distances, "Inter-Host Distances")

    # Statistical tests
def perform_statistical_test(group1, group2, group1_label, group2_label):
   group1_normal = lilliefors(group1)[1] > 0.05 
   group2_normal = lilliefors(group2)[1] > 0.05
   if group1_normal and group2_normal:
       # T-test if both groups are normally distributed
       t_stat, p_value = stats.ttest_ind(group1, group2, equal_var=False)
       test_name = 'T-test'
   else:
       # Mann-Whitney U test if at least one group is not normally distributed
       u_stat, p_value = stats.mannwhitneyu(group1, group2, alternative='two-sided')
       test_name = 'Mann-Whitney U test'

       #Determine Significance
       significance = 'Significant' if p_value < 0.05 else 'Non-Significant'
   print(f"{test_name} p-value ({group1_label} vs {group2_label}): {p_value} ({significance})")
print("\nStatistical Tests for Inter-host Distances:")

# Compare inter-host between Bovine and Human  
for host, distances in intra_distances.items():
    perform_statistical_test(intra_distances['Bovine'], distances, "Bovine Intra-host Distance",  f"{host} Intra-host Distances") 
    perform_statistical_test(intra_distances['Human'], distances, "Human Intra-host Distance",  f"{host} Intra-host Distances") 
    perform_statistical_test(inter_distances, distances, "Inter-host Distance", f"{host} Intra-host Distances")
#Output Intra-Host SNP Variation:
print("\nAverage Intra-host SNP Distance:")
for host, avg_distance in avg_intra_distances.items():
   print(f"{host}: {avg_distance:.2f}")

#Output Inter-Host SNP Variation:
print("\nAverage Inter-host SNP Distance:")
print(f"{avg_inter_distance:.2f}")

#Determine if the groups are normal can explain if T test or Mann-Whitney was used for analysis
print ("\nNormality Results:")
for host, is_normal in normality_results.items():
    print(f"{host} {'normally distributed' if is_normal else 'not normally distributed'}")


