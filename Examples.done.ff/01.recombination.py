"""
Original Author: Yuwen Pan
Modified by: Xiaoxi Zhang (To adapt this script for the pipeline)

This script calculates the recombination rate based on genomic marker positions and outputs the results to a specified file.

Input:
    1. A tab-separated file with no header. Columns are assumed to represent:
        - Column 1: Start position of genomic markers
        - Column 2: Genomic marker identifier
        - Column 3: Genetic map distance (in cM)
Output:
    - A space-separated file containing:
        - Start position
        - Recombination rate per base pair
        - Recombination rate in cM per megabase
        - Genetic map distance
"""

import pandas
import sys
import numpy as np

# Retrieve input and output file paths from command-line arguments
infile = sys.argv[1]
recfile = sys.argv[2]

# Read the input file, assuming it is tab-separated with no header
inf = pandas.read_csv(infile, sep='\t', header=None)

# Set the DataFrame index to the second column (genomic marker identifier)
inf.index = inf[1]

# Create a new DataFrame to store recombination rate data
rec = pandas.DataFrame()
rec['start.pos'] = list(inf[1])
rec['recom.rate.perbp'] = 0.0
rec['recom.rate.cMperMb'] = 0.0
rec['gm'] = list(inf[2])
rec['next_gm'] = rec['gm'].shift(periods=-1)

# Set the last row of 'next_gm' to its own 'gm' value to avoid NaN
rec.iloc[-1, -1] = rec.iloc[-1, -2]
rec['next_pos'] = rec['start.pos'].shift(periods=-1)

# Set the last row of 'next_pos' to 0 to avoid division by zero
rec.iloc[-1, -1] = 0.0

# Calculate recombination rates
rec['recom.rate.perbp'] = ((rec['next_gm'] - rec['gm']) / (rec['next_pos'] - rec['start.pos'])) / 100.0       #Mperbp
rec['recom.rate.cMperMb'] = (rec['next_gm'] - rec['gm']) / ((rec['next_pos'] - rec['start.pos']) / 1000000 )  #cMperMb

# Drop unnecessary columns
rec.drop(['next_gm', 'next_pos'], axis=1, inplace=True)

# Set the last row of 'recom.rate.perbp' and 'recom.rate.cMperMb' to 0
rec.iloc[-1, -2] = 0.0
rec.iloc[-1, -3] = 0.0

# Write the results to the output file
rec.to_csv(recfile, sep=' ', index=None, float_format='%.15f')

