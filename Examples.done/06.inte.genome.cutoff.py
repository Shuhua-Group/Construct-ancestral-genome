"""
This script integrates ancestral gene pool across multiple files to generate an ancestral haplotype.

Author: Xiaoxi Zhang
Purpose: Integrate ancestral gene pools and filter variants based on user-defined cutoffs to generate an ancestral haplotype.

Input:
    1. Output ancestral haplotype (gzipped).
    2. Missing rate cutoff: Minimum missing rate (e.g., 0.2).
    3. Length cutoff: The maximum length that can be extended at one time when obtaining ancestral segments on a haplotype (e.g., 10000.0 bp).
    4. Input files: One or more gzipped files containing ancestral gene pool.

Output:
    - A gzipped file containing filtered genomic data based on the specified cutoffs.
"""

import sys,os,gzip,random
# Input parameters
outfile = sys.argv[1]           # Output file path
cutoff_freq = float(sys.argv[2])  # Missing rate cutoff (e.g., 0.2)
cutoff_length = float(sys.argv[3])  # Length cutoff (e.g., 10000.0)
infiles = sys.argv[4:]          # Input gzipped ancestral gene pool files

# Load and integrate ancestral gene pool from input files
info = []  # List to store integrated ancestral gene pool
for infile in infiles:
	k = 0
	f = gzip.open( infile, 'rt' )
	f.readline()   # Skip the header line
	for line in f:
		y = line.split()
		chr, pos, genos = y[0], y[1], y[2:]   # Chromosome, position, genotypes
		if k + 1 > len(info):
			info.append( [ chr, pos, '|'.join(genos) ] )
			k += 1
		else:
			info[k][2] = info[k][2] + '|' + '|'.join(genos)   # Append genotypes
			k += 1
	f.close()

# Function to handle random selection for missing data
def deal_one_row(point, n, slgeno):
	"""
	Randomly select a valid allele from the genotype list.
	If all genotypes are missing, return '.' and a random index.
	"""
	d = { point: '' }
	while len(d) < n:
		t = random.randint( 0, n-1 )
		if t not in d:
			if slgeno[t] != '.':
				return slgeno[t], t
			else:
				d[t] = ''
		else:
			continue
	return '.', random.randint( 0, n-1 )

# Function to judge if a row passes the frequency cutoff
def judge( li ):
	"""
	Check if the frequency of non-missing genotypes meets the cutoff.
	"""
	kl = 0
	for i in range( len(li) ):
		if li[i] != '.':
			kl += 1
	if cutoff_freq <= float(kl) / len(li):
		return 'ok'
	else:
		return 'no'

# Write filtered data to output file
fout = gzip.open( outfile, 'wt' )
n = len(info[0][2].split('|'))   # Number of samples
point = random.randint( 0, n-1 )    # Random starting point for selection
start = '.'
for i in range( len(info) ):
	slgeno = info[i][2].split('|')   # Split genotypes for the current row
	res = judge( slgeno )   # Check frequency cutoff
	if res != 'ok':   # If frequency cutoff is not met
		point = random.randint( 0, n-1 )
		fout.write( info[i][0] + '\t' + info[i][1] + '\t' + '.' + '\n' )
		start = '.'
		continue
	if slgeno[point] == '.':   # Handle missing genotype at the current point
		fw, point = deal_one_row( point, n, slgeno )
		if fw == '.':
			fout.write( info[i][0] + '\t' + info[i][1] + '\t' + '.' + '\n' )
			start = '.'
		else:
			fout.write( info[i][0] + '\t' + info[i][1] + '\t' + fw  + '\n' )
			start = info[i][1]
	else:   # If the genotype is valid
		result = ''
		if start == '.':
			result = 'ok'
		else:
			if int(info[i][1]) + 1 - int(start) > cutoff_length:   # Check length cutoff
				result = 'no'
		if result == 'no':
			fw, point = deal_one_row( point, n, slgeno )
			if fw == '.':
				fout.write( info[i][0] + '\t' + info[i][1] + '\t' + fw  + '\n' )
				start = '.'
			else:
				fout.write( info[i][0] + '\t' + info[i][1] + '\t' + fw  + '\n' )
				start = info[i][1]
		else:
			fout.write( info[i][0] + '\t' + info[i][1] + '\t' + slgeno[point] + '\n' )
fout.close()
