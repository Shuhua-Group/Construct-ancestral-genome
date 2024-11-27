"""
This script extracts identified ancestral alleles. 
It extracts alleles belonging to the specified ancestry from the selected samples in phased VCF files and outputs an ancestral gene pool file.

Author: Xiaoxi Zhang
Purpose: To extract and organize identified ancestral alleles for downstream analysis.

Input:
    1. Input file (gzipped): Contains ancestry information data to process.
    2. Sample list: A file with sample IDs to select for analysis.
    3. Phased VCF file (gzipped): Contains variant data with ancestral allele information.
    4. Output file (gzipped): Output an ancestral gene pool file.
    5. Sources: Additional sources of ancestral allele information.

Output:
    - A gzipped file with extracted ancestral allele information for the selected samples.
"""

import sys,os,gzip

# Input and output files
infile = sys.argv[1]  # Gzipped input file with ancestry information data
inlist = sys.argv[2]  # Sample list file
invcf = sys.argv[3]   # Gzipped VCF file with variant data
outfile = sys.argv[4] # Gzipped output file for processed results
sources = sys.argv[5:] # Additional sources of ancestral allele information

# Load sample list
d = {}   # Dictionary to store selected sample IDs
f = open( inlist )
for line in f:
	y = line.split()
	d[ y[0] ] = ''
f.close()

# Load additional sources
sourc = {}
for source in sources:
	sourc[ source ] = ''

# Parse the VCF file to extract ancestral allele information
info = {}       # Dictionary to store ancestral allele data
id_rank = {}    # Mapping between sample ID and haplotypes
select_id = {}  # Selected sample indices
f = gzip.open( invcf, 'rt' )
for line in f:
	if line[0] == '#':
		y = line.split()
		if y[0] == '#CHROM':   # Parse header line to identify samples
			k = 0
			for i in range( 9, len( y ) ):
				if y[i] in d:
					select_id[ i ] = ''
					id_rank[ y[i] + '_1' ] = k
					id_rank[ y[i] + '_2' ] = k + 1
					k += 2
			continue
		else:
			continue
	tttt = []
	y = line.split()
	chr, pos = y[0], y[1]
	if chr[:3] == 'chr':
		chr = chr[3:]   # Remove "chr" prefix if present
	for i in range( 9, len(y) ):
		if i in select_id:
			tttt.append( y[i][0] )   # Extract allele for haplotype 1
			tttt.append( y[i][2] )   # Extract allele for haplotype 2
	info[ chr + ':' + pos ] = ''.join( tttt )
f.close()

# Handle the header (title) line of the input file
def deal_title( t ):
	"""
	Parse the title row to extract sample information and create a mapping.
	"""
	fr = {}
	chr, pos, ids = t[0], t[1], t[2:]
	for i in range( 0, len( ids ) ):
		tru_id, tru_hap = ids[i].split('_')[0], ids[i].split('_')[1]	
		if tru_id in d:
			fr[i+2] = [ tru_id, tru_hap ]
		else:
			print( 'Not.include' )
	return fr

# Process a single row of data
def deal_one_row( y ):
	"""
	Process a single row of input file to extract ancestral allele information.
	"""
	chr, pos, results = y[0], y[1], y[2:]
	tem = [  ]
	for i in range( 2, len(y) ):
		if i in re:
			if y[i] in sourc:
				i_id, i_hap = re[i][0], re[i][1]
				allele = info[ chr + ':' + pos ][ id_rank[ i_id + '_' + i_hap ] ]
				tem.append( allele )
			else:
				tem.append( '.' )
	return tem

# Process input file and write output
fout = gzip.open( outfile, 'wt' )
f = gzip.open( infile, 'rt' )
title = f.readline().split()
re = deal_title( title ) # Create mapping for samples
# Write header to output file
ti = [ 'Chr', 'Pos' ]
for r in re:
	ti.append( '_'.join(re[r]) )
fout.write( '\t'.join( ti ) + '\n' )
# Write rows of processed data
for line in f:
	y = line.split()
	chr, pos, results = y[0], y[1], y[2:]
	fw = deal_one_row( y )
	fout.write( '\t'.join( [ chr, pos ] + fw ) + '\n' )
f.close()
fout.close()
