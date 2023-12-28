#!/usr/bin/env python3

# Released under the MIT license.

# Split regions in BED into n files with approximately equal region sizes. 
# A region is never split. 13 is a good number. 

import sys
import pandas as pd

chromosome_data = pd.read_csv(sys.argv[1], names = ['chr', 'start', 'stop'], usecols=range(3), sep = '\t')

chromosome_data['size'] = chromosome_data['stop'] - chromosome_data['start']

# Number of bins
n = int(sys.argv[2])

# Sort chromosome data by size in descending order
sorted_data = chromosome_data.sort_values(by='size', ascending=False)

# Initialize empty bins as lists
bins = [[] for _ in range(n)]

# Allocate chromosomes to bins
for index, row in sorted_data.iterrows():
    # Find the bin with the fewest chromosomes
    min_bin = min(range(n), key=lambda i: sum(chrom['size'] for chrom in bins[i]))

    # Place the chromosome data in the selected bin
    bins[min_bin].append(row.to_dict())

# Create a DataFrame to store the results
result_df = pd.DataFrame({
    'bin': [i + 1 for i in range(n) for _ in bins[i]],
    'chr': [chromosome['chr'] for bin_chromosomes in bins for chromosome in bin_chromosomes],
    'start': [int(chromosome['start']) for bin_chromosomes in bins for chromosome in bin_chromosomes],
    'stop': [int(chromosome['stop']) for bin_chromosomes in bins for chromosome in bin_chromosomes],
    'size': [int(chromosome['size']) for bin_chromosomes in bins for chromosome in bin_chromosomes]
})

# Print the result DataFrame, ordered by size within each bin
result_df = result_df.sort_values(by=['bin', 'size'], ascending=[True, False])

for id, group in result_df.groupby(['bin']):
    group[['chr', 'start', 'stop']].to_csv(f'{id}.bed', index=False, header=False, sep = '\t')
