"""
This script identifies the ancestry of haplotypes based on ChromoPainter output.
It processes a list of haplotype segments and assigns ancestry labels based on the majority vote for each segment.

Author:  Xiaoxi Zhang
Purpose: To adapt ChromoPainter output for downstream reconstruction of ancestral genomes.

Input:
    1. ChromoPainter output file (gzipped). Each haplotype is assumed to be painted 10 times by default in the ChromoPainter output.
    2. Population information file (tab-separated):
        - Column 1: Individual ID
        - Column 2: Population label
        - Column 3: Rank or index used in ChromoPainter.
    3. Output file path (gzipped).
    4. Cutoff value: Minimum number of votes required to assign ancestry.

Output:
    - A gzipped file containing haplotypes annotated with ancestry labels or '.' if no label meets the cutoff.
"""


import sys,os,gzip

# Input file paths and cutoff value
infile  = sys.argv[1]   # ChromoPainter output file
inlist  = sys.argv[2]   # ID, population, rank information
outfile = sys.argv[3]   # Output file path
cutoff  = int(sys.argv[4])   # Minimum votes required for ancestry assignment (>=cutoff)


# Load population and rank mapping
d = {}
f = open( inlist )
for line in f:
	y = line.split()
	id, pop, rank = y[0], y[1], y[2]
	d[ rank ] = [ id, pop ]
f.close()


# Process a single row of ChromoPainter data
def deal_one_row( se, y, inf ):
	"""
	Processes one row of ChromoPainter output.
	Accumulates population counts for each haplotype position.
	"""
	for k in range( 1, len(y) ):
		pop = d[ y[k] ][1]   # Population label
		if k not in se:
			se[ k ] = { pop: 1 }   # Initialize population count
		else:
			if pop not in se[k]:
				se[ k ][ pop ] = 1   # Add new population
			else:
				se[ k ][ pop ] += 1  # Increment population count
	return se

# Process a single haplotype
def deal_one_hap( u ):
	"""
	Processes one haplotype (assumes 11 rows per haplotype).
	Determines ancestry label for each position based on majority vote.
	"""
	inf = u[0].split()   # Sample information row
	se = {}
	for j in range( 1, len(u) ):
		y = u[j].split()
		se = deal_one_row( se, y, inf )

	for mem in se:
		# Sort populations by count in descending order
		p_sort = sorted( se[mem].items(), key = lambda x:x[1], reverse = True )
		if p_sort[0][1] >= cutoff:
			inf.append( p_sort[0][0] )   # Assign ancestry label
		else:
			inf.append( '.' )   # No label meets cutoff
	return inf
	

# Process the input file and write the results
fout = gzip.open( outfile, 'wt' )
l = []
f = gzip.open( infile, 'rt' )
f.readline()# Skip the header line
lines = f.readlines()
f.close()
n = len(lines)
for i in range( 0, int(n/11) ):   # Process haplotypes in blocks of 11 lines
	start, end = 11 * i, 11 * i + 10
	fw = deal_one_hap( lines[ start: end + 1 ] )   # Process 11 lines per haplotype
	fout.write( '\t'.join( fw ) + '\n' )   # Write output
fout.close()
