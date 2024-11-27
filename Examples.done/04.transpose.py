"""
This script transposes ancestry information from a positional format to a sample-centric format.
It combines positional information with sample-specific data to produce a table with samples as columns.

Author: Xiaoxi Zhang
Purpose: Transform and reformat ancestry information data for downstream analysis.

Input:
    1. A gzipped input file containing ancestry information.
    2. A position list file (tab-separated):
        - Column 1: Chromosome
        - Column 2: Position
    3. Output file path and name (gzipped).

Output:
    - A gzipped output file with the following format:
        - Header: Chr, Pos, followed by sample identifiers (Sample_Pos format).
        - Data: Transposed positional data for each sample.
"""

import sys,os,gzip

# Input file paths
infile  = sys.argv[1]   # Gzipped input file with ancestry information
poslist = sys.argv[2]   # Position list file
outfile = sys.argv[3]   # Output file path and name (gzipped)

# Read position list file
d = {}   # Dictionary to store chromosome and position information
k = 1   # Counter for position index
f = open( poslist )
for line in f:
	y = line.split()
	chr, pos = y[0], y[2]   # Extract chromosome and position
	d[ k ] = [ chr, pos ]   # Store in dictionary
	k += 1
f.close()

# Read ancestry information data from gzipped input file
l = []   # List to store rows of the input file
f = gzip.open( infile, 'rt' )
for line in f:
	y = line.split()   # Split each line into columns
	l.append( y )   # Append to list
f.close()

# Write output file
fout = gzip.open( outfile, 'wt' )
# Write header
fout.write( 'Chr\tPos' )   # Fixed columns for chromosome and position
for i in range( len(l) ):
	fout.write( '\t' + l[i][2] + '_' + l[i][1] )   # Sample identifier (Sample_Pos format)
fout.write( '\n' )

# Write transposed data
for i in range(1,len(d)+1):   # Loop through positions
	t = d[i]   # Retrieve chromosome and position
	for j in range( len(l) ):   # Loop through samples
		t.append( l[j][2+i] )   # Append sample-specific data for the position
	fout.write( '\t'.join( t ) + '\n' )   # Write row to output
fout.close()
