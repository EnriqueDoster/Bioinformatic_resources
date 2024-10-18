import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

# Configuration
output_suffix = '0.0005'

base_input_directory = f"Temp_dir_{output_suffix}"
metadata_filename = f"output_TreeCluster_labels_avg_clade_t_{output_suffix}.txt"  # mh_genome_groups.txt
output_file_path = f"merged_kraken_output_results_{output_suffix}.txt"
output_sn_sp_path = f"sn_sp_precision_by_confidence_{output_suffix}.txt"
output_figure = f"figure_sn_sp_precision_by_confidence_{output_suffix}.png"  # Assuming it's a PNG file, not txt


#base_input_directory = 'Temp_dir'
#metadata_filename = 'mh_genome_groups.txt' # .txt
#output_file_path = 'merged_output_results.txt'


# Helper function to process filenames
def process_filename(filename):
    return filename.replace('modified_', '').split('.')[0]

# Read metadata and create mappings, adding "NonMh" for missing entries
metadata_df = pd.read_csv(metadata_filename, delimiter='\t', header=None, names=['SequenceName', 'ProposedLabel'])
metadata_df['SequenceName'] = metadata_df['SequenceName'].apply(process_filename)
metadata_map = dict(zip(metadata_df['SequenceName'], metadata_df['ProposedLabel']))

# Initialize a list to hold all concatenated data
all_data = []

# Iterate through each "Run_" directory
for run_directory in glob.glob(os.path.join(base_input_directory, "Run_*")):
    run_name = os.path.basename(run_directory)
    report_files_pattern = os.path.join(run_directory, "*.report*")

    for report_file in glob.glob(report_files_pattern):
        # Read report file with specified column names
        df = pd.read_csv(report_file, delimiter='\t', header=None, 
                         names=["Perc_frag_rooted", "Num_frag_rooted", "Num_frag_taxon", 
                                "Taxa_rank", "NCBI_ID", "Taxa_name"])
        original_filename = os.path.basename(report_file).replace('.report', '')
        processed_filename = process_filename(original_filename)
        # Add extra columns for processed filename, original filename, and run name
        df['SequenceName'] = processed_filename
        df['FileName'] = original_filename
        df['Run'] = run_name
        # Map to ProposedLabel and fill missing entries with "NonMh"
        df['ProposedLabel'] = df['SequenceName'].map(metadata_map).fillna("NonMh")
       #print(df)
        all_data.append(df)

# Concatenate all data from different runs into a single DataFrame
concatenated_data = pd.concat(all_data, ignore_index=True)

# Save concatenated data to file
concatenated_data.to_csv(output_file_path, sep='\t', index=False)

##
### This section is to calculate whether the correct taxa was identified for each sample
##

# Analyze and output the final results
unique_files = concatenated_data['FileName'].unique()
results = []

# Loop through the results for each unique file
for file in unique_files:
    file_data = concatenated_data[concatenated_data['FileName'] == file]
    group = file_data['ProposedLabel'].iloc[0]

    # Extract all rows matching "Mh_" pattern
    mh_rows = file_data[file_data['Taxa_name'].str.contains(r"\s*Mh_\S+", regex=True)]

    # Check for an exact match with "Mh_{group}"
    exact_matches = mh_rows[mh_rows['Taxa_name'].str.contains(fr"\s*Mh_{group}\b", regex=True)]
    match_found = 'Y' if not exact_matches.empty else 'N'

    # OnTargetReads: Sum of 'Num_frag_taxon' for exact matches
    OnTargetReads = exact_matches['Num_frag_taxon'].sum() if not exact_matches.empty else 0

    # TotalReads: Sum of 'Num_frag_taxon' for all rows in file_data
    TotalReads = file_data['Num_frag_taxon'].sum()

    # Check for multiple "Mh_" matches and calculate OffTargetPSVReads
    multiple_matches = 'Y' if len(mh_rows) > 1 else 'N'
    # OffTargetPSVReads: Sum of 'Num_frag_taxon' for all mh_rows except for the exact matches
    OffTargetPSVReads = mh_rows['Num_frag_taxon'].sum() - OnTargetReads if multiple_matches == 'Y' else 0

    # Extract "Confidence" and "Run" from the filename
    confidence_match = re.search(r'conf_(\d+\.\d+)', file)  # Adjusted regex to be more generic
    run_match = re.search(r'run_(\d+)', file)

    confidence = confidence_match.group(1) if confidence_match else "Unknown"
    run_number = run_match.group(1) if run_match else "Unknown"

    results.append([file, group, match_found, multiple_matches, confidence, run_number, OnTargetReads, TotalReads, OffTargetPSVReads])

# Convert the results list to a DataFrame and adjust column names accordingly
results_df = pd.DataFrame(results, columns=['SequenceName', 'Group', 'MatchFound', 'MultipleMatches', 'Confidence', 'Run', 'OnTargetReads', 'TotalReads', 'OffTargetPSVReads'])

results_df.to_csv(output_file_path, sep='\t', index=False)
#print(results_df)

print(f"Results saved to {output_file_path}")

##
### The next part is optional to graph the results
##

# Initialize a list to store the results summary
results_summary = []

# Iterate over each unique combination of Run and Confidence
for run in sorted(results_df['Run'].unique()):
    for confidence in sorted(results_df['Confidence'].unique()):
        subset = results_df[(results_df['Run'] == run) & (results_df['Confidence'] == confidence)]
        # Calculate True Positives (TP), False Negatives (FN), False Positives (FP), and True Negatives (TN)
        TP = len(subset[(subset['Group'] != 'NonMh') & (subset['MatchFound'] == 'Y')])
        FN = len(subset[(subset['Group'] != 'NonMh') & (subset['MatchFound'] == 'N')])
        FP = len(subset[(subset['Group'] == 'NonMh') & ((subset['MultipleMatches'] == 'Y') | (subset['MatchFound'] == 'Y'))])
        TN = len(subset[(subset['Group'] == 'NonMh') & ((subset['MultipleMatches'] == 'N') & (subset['MatchFound'] == 'N'))])

        # Adjust TP for Precision calculation
        Precision_TP = len(subset[(subset['Group'] != 'NonMh') & (subset['MatchFound'] == 'Y')  & (subset['MultipleMatches'] == 'N')])
        Total_PSV = len(subset[(subset['Group'] != 'NonMh')])
                
        # Calculate Sensitivity and Specificity
        sensitivity = TP / (TP + FN) if (TP + FN) > 0 else 0
        specificity = TN / (TN + FP) if (TN + FP) > 0 else 0
        
        # Calculate Precision
        precision = Precision_TP / Total_PSV if (Total_PSV) > 0 else 0

        results_summary.append([run, confidence, sensitivity, specificity,precision])

summary_df = pd.DataFrame(results_summary, columns=['Run', 'Confidence', 'Sensitivity', 'Specificity','Precision'])


## Optionally, save summary_df to a CSV file
summary_df.to_csv(output_sn_sp_path, index=False)
print(f"Saved summary to {output_sn_sp_path}")

# Calculate the average and standard deviation for sensitivity, specificity, and precision at each confidence level across all runs
avg_std_metrics = summary_df.groupby('Confidence').agg({
    'Sensitivity': ['mean', 'std'],
    'Specificity': ['mean', 'std'],
    'Precision': ['mean', 'std']
}).reset_index()

# Flatten the MultiIndex for columns created by aggregation
avg_std_metrics.columns = ['Confidence', 'Sensitivity_mean', 'Sensitivity_std', 
                           'Specificity_mean', 'Specificity_std', 'Precision_mean', 'Precision_std']

# Plotting
plt.figure(figsize=(12, 6))

# Sensitivity
plt.errorbar(avg_std_metrics['Confidence'], avg_std_metrics['Sensitivity_mean'], 
             yerr=avg_std_metrics['Sensitivity_std'], fmt='-o', capsize=5, 
             label='Average Sensitivity', color='blue')

# Specificity
plt.errorbar(avg_std_metrics['Confidence'], avg_std_metrics['Specificity_mean'], 
             yerr=avg_std_metrics['Specificity_std'], fmt='-x', capsize=5, 
             label='Average Specificity', color='red')

# Precision
plt.errorbar(avg_std_metrics['Confidence'], avg_std_metrics['Precision_mean'], 
             yerr=avg_std_metrics['Precision_std'], fmt='-s', capsize=5, 
             label='Average Precision', color='green')

plt.title('Average Sensitivity, Specificity, and Precision by Confidence Level')
plt.xlabel('Confidence Level')
plt.xticks(rotation=45)  # Rotate x-axis labels for better readability
plt.ylabel('Metric Value')
plt.legend()
plt.tight_layout()  # Adjust layout to make room for the rotated x-axis labels
plt.grid(True)

# Save the plot
plt.savefig(output_figure)
print(f"Plot saved as {output_figure}.")
