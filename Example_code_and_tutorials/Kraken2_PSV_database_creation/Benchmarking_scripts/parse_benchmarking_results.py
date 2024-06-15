import pandas as pd
import glob
import os
import re
import matplotlib.pyplot as plt

# This script works on a single output run with multiple kraken report files

# Input directory and metadata file path
input_directory = 'Temp_dir/Run_1'
metadata_filename = 'mh_genome_groups.txt'

# Output file path
output_file_path = 'merged_output_results.txt'

# Pattern to find report files
report_files_pattern = os.path.join(input_directory, "*.report*")

# Function to process filenames
def process_filename(filename):
    # Remove 'modified_' prefix and everything after the first '.'
    processed_name = filename.replace('modified_', '').split('.')[0]
    return processed_name

# Read metadata, process 'SequenceName', and create a mapping
metadata_df = pd.read_csv(metadata_filename, delimiter='\t', header=None, names=['SequenceName', 'ProposedLabel'])
metadata_df['SequenceName'] = metadata_df['SequenceName'].apply(process_filename)
metadata_map = dict(zip(metadata_df['SequenceName'], metadata_df['ProposedLabel']))

# Process each report file
all_data = []


for report_file in glob.glob(report_files_pattern):
    df = pd.read_csv(report_file, delimiter='\t', header=None)
    # Extract the original basename of the file (without the path)
    original_filename = os.path.basename(report_file).replace('.report', '')
    # Process the report filename for matching with metadata
    processed_filename = process_filename(original_filename)
    # Insert both processed and original filenames as the first and second columns
    df.insert(0, 'SequenceName', processed_filename)
    df.insert(1, 'FileName', original_filename)
    all_data.append(df)

# Concatenate all data
concatenated_data = pd.concat(all_data, ignore_index=True)

# Map processed SequenceName to ProposedLabel
concatenated_data['ProposedLabel'] = concatenated_data['SequenceName'].map(metadata_map)
#print(concatenated_data)

#concatenated_data_df = pd.DataFrame(concatenated_data)
#concatenated_data_df.to_csv('temp_concat_data.tsv', sep='\t', index=False)

# Analyze and output the final results
unique_files = concatenated_data['FileName'].unique()
results = []

# Loop though the results for each unique file
for file in unique_files:
    file_data = concatenated_data[concatenated_data['FileName'] == file]
    group = file_data['ProposedLabel'].iloc[0]
    
    # Extract all rows matching "Mh_" pattern
    mh_rows = file_data[file_data.iloc[:, -2].str.contains(r"\s*Mh_\S+", regex=True)]
    
    # Check for an exact match with "Mh_{group}"
    exact_matches = mh_rows[mh_rows.iloc[:, -2].str.contains(fr"\s*Mh_{group}\b", regex=True)]
    match_found = 'Y' if not exact_matches.empty else 'N'
    
    # Check for multiple "Mh_" matches
    multiple_matches = 'Y' if len(mh_rows) > len(exact_matches) else 'N'
    
    # Extract "Confidence" from the filename
    # Assuming filename format: something_conf_X.Y.txt where X.Y is the confidence level
    start = file.find('conf_') + 5  # +5 to skip past 'conf_'
    end = file.find('.txt', start)
    confidence = file[start:end] if start > 4 else "Unknown"  # Check if 'conf_' was found
    
    results.append([file, group, match_found, multiple_matches, confidence])

# Convert the results list to a DataFrame and save it to a text file
results_df = pd.DataFrame(results, columns=['SequenceName', 'Group', 'MatchFound', 'MultipleMatches', 'Confidence'])
results_df.to_csv(output_file_path, sep='\t', index=False)

print(f'Results saved to {output_file_path}')


##
### The next part is optional to graph the results
##

# Assuming results_df is already defined
results_summary = []

for confidence in sorted(results_df['Confidence'].unique()):
    subset = results_df[results_df['Confidence'] == confidence]
    
    # Sensitivity: Proportion of actual positives that were correctly identified
    match_found_count = len(subset[subset['MatchFound'] == 'Y'])
    total_possible_matches = len(subset)  # Total rows that could have a match
    sensitivity = match_found_count / total_possible_matches if total_possible_matches > 0 else 0
    
    # Precision: Proportion of positive identifications that were actually correct
    correct_and_unique_matches = len(subset[(subset['MatchFound'] == 'Y') & (subset['MultipleMatches'] == 'N')])
    precision = correct_and_unique_matches / match_found_count if match_found_count > 0 else 0
    
    results_summary.append([confidence, precision, sensitivity])

# Convert summary to DataFrame
summary_df = pd.DataFrame(results_summary, columns=['Confidence', 'Precision', 'Sensitivity'])

# Plotting Precision and Sensitivity vs. Confidence
plt.figure(figsize=(10, 6))
plt.plot(summary_df['Confidence'], summary_df['Precision'], label='Precision', marker='o', linestyle='-', color='blue')
plt.plot(summary_df['Confidence'], summary_df['Sensitivity'], label='Sensitivity', marker='x', linestyle='--', color='red')

plt.title('Precision and Sensitivity by Confidence Level')
plt.xlabel('Confidence Level')
plt.ylabel('Metric Value')
plt.legend()
plt.grid(True)

# Save the plot
plt.savefig('precision_sensitivity_by_confidence.png')
print("Plot saved as precision_sensitivity_by_confidence.png")