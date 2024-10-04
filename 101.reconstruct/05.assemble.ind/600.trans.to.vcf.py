import sys,os,gzip

invcf  = sys.argv[1]
intxt1 = sys.argv[2]
intxt2 = sys.argv[3]
samID  = sys.argv[4]
outfil = sys.argv[5]

d = {}
f = gzip.open( invcf, 'rt' )
for line in f:
	if line[0] == '#':
		continue
	y = line.split()
	chr, pos, ref, alt = y[0], y[1], y[3], y[4]
	if chr[:3] == 'chr':
		chr = chr[3:]
	d[ chr + ':' + pos ] = [ ref, alt ]
f.close()


f = gzip.open( intxt1, 'rt' )
for line in f:
	y = line.split()
	chr, pos, allele = y[0], y[1], y[2][0]
	d[ chr + ':' + pos ].append( allele )
f.close()


f = gzip.open( intxt2, 'rt' )
for line in f:
	y = line.split()
	chr, pos, allele = y[0], y[1], y[2][0]
	d[ chr + ':' + pos ].append( allele )
f.close()


fout = open( outfil, 'w' )
fout.write( '##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description="All filters passed">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##contig=<ID=1>\n##contig=<ID=2>\n##contig=<ID=3>\n##contig=<ID=4>\n##contig=<ID=5>\n##contig=<ID=6>\n##contig=<ID=7>\n##contig=<ID=8>\n##contig=<ID=9>\n##contig=<ID=10>\n##contig=<ID=11>\n##contig=<ID=12>\n##contig=<ID=13>\n##contig=<ID=14>\n##contig=<ID=15>\n##contig=<ID=16>\n##contig=<ID=17>\n##contig=<ID=18>\n##contig=<ID=19>\n##contig=<ID=20>\n##contig=<ID=21>\n##contig=<ID=22>\n##contig=<ID=X>\n##contig=<ID=Y>\n##contig=<ID=MT>\n##contig=<ID=M>\n' )
fout.write( '\t'.join( '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT'.split() + [ samID ] ) + '\n' )
for dd in d:
	chr, pos, ref, alt, a1, a2 = dd.split(':')[0], dd.split(':')[1], d[dd][0], d[dd][1], d[dd][2], d[dd][3]
	fout.write( '\t'.join( [ chr, pos, '.', ref, alt, '.', 'PASS', '.', 'GT', a1 + '|' + a2 ] ) + '\n' )
fout.close()

os.system( 'rm ' + outfil + '.gz' )
os.system( 'bgzip -f ' + outfil )
os.system( 'tabix -f -p vcf ' + outfil + '.gz' )
