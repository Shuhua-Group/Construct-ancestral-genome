"""
This script converts two ancestral haplotypes into a VCF-format ancestral individual. 
It integrates reference and alternate allele information from input files and formats it for downstream analysis.

Author: Xiaoxi Zhang
Purpose: Generate an ancestral individual using two ancestral haplotypes.

Input:
    1. invcf: A gzipped VCF file with reference and alternate allele information.
    2. intxt1: A gzipped file containing the first ancestral haplotype.
    3. intxt2: A gzipped file containing the second ancestral haplotype.
    4. samID: Specify the name of the generated ancestral sample.
    5. outfil: Output VCF file path.

Output:
    - A compressed VCF file with genotype information.
    - Associated index file for the VCF file.
"""

import sys, os, gzip

# Input and output file paths
invcf = sys.argv[1]   # Input gzipped VCF file
intxt1 = sys.argv[2]  # Ancestral haplotype 1
intxt2 = sys.argv[3]  # Ancestral haplotype 2
samID = sys.argv[4]   # Sample ID
outfil = sys.argv[5]  # Output VCF file

# Step 1: Parse the reference and alternate alleles from the input VCF file
d = {}  # Dictionary to store reference and alternate allele information
f = gzip.open( invcf, 'rt' )
for line in f:
	if line[0] == '#':   # Skip header line
		continue
	y = line.split()
	chr, pos, ref, alt = y[0], y[1], y[3], y[4]
	if chr[:3] == 'chr':   # Remove "chr" prefix if present
		chr = chr[3:]
	d[ chr + ':' + pos ] = [ ref, alt ]
f.close()

# Step 2: Add allele information from intxt1
f = gzip.open( intxt1, 'rt' )
for line in f:
	y = line.split()
	chr, pos, allele = y[0], y[1], y[2][0]   # Extract chromosome, position, and allele
	d[ chr + ':' + pos ].append( allele )
f.close()

# Step 3: Add allele information from intxt2
f = gzip.open( intxt2, 'rt' )
for line in f:
	y = line.split()
	chr, pos, allele = y[0], y[1], y[2][0]   # Extract chromosome, position, and allele
	d[ chr + ':' + pos ].append( allele )
f.close()

# Step 4: Write VCF output
fout = open( outfil, 'w' )
# Write VCF header
fout.write( '##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description="All filters passed">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##contig=<ID=1>\n##contig=<ID=2>\n##contig=<ID=3>\n##contig=<ID=4>\n##contig=<ID=5>\n##contig=<ID=6>\n##contig=<ID=7>\n##contig=<ID=8>\n##contig=<ID=9>\n##contig=<ID=10>\n##contig=<ID=11>\n##contig=<ID=12>\n##contig=<ID=13>\n##contig=<ID=14>\n##contig=<ID=15>\n##contig=<ID=16>\n##contig=<ID=17>\n##contig=<ID=18>\n##contig=<ID=19>\n##contig=<ID=20>\n##contig=<ID=21>\n##contig=<ID=22>\n##contig=<ID=X>\n##contig=<ID=Y>\n##contig=<ID=MT>\n##contig=<ID=M>\n' )
fout.write( '\t'.join( '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT'.split() + [ samID ] ) + '\n' )
# Write variant information
for dd in d:
	chr, pos, ref, alt, a1, a2 = dd.split(':')[0], dd.split(':')[1], d[dd][0], d[dd][1], d[dd][2], d[dd][3]   # Extract chromosome, position and  alleles
	fout.write( '\t'.join( [ chr, pos, '.', ref, alt, '.', 'PASS', '.', 'GT', a1 + '|' + a2 ] ) + '\n' )
fout.close()

# Step 5: Compress and index the VCF file
os.system( 'bgzip -f ' + outfil )   # Compress the VCF file
os.system( 'tabix -f -p vcf ' + outfil + '.gz' )   # Create an index for the VCF file
