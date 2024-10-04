import sys,os,gzip

infile  = sys.argv[1]
inlist  = sys.argv[2]
invcf   = sys.argv[3]
outfile = sys.argv[4]
sources = sys.argv[5:]


d = {}
f = open( inlist )
for line in f:
	y = line.split()
	d[ y[0] ] = ''
f.close()

sourc = {}
for source in sources:
	sourc[ source ] = ''


info = {}
id_rank = {}
select_id = {}
f = gzip.open( invcf, 'rt' )
for line in f:
	if line[0] == '#':
		y = line.split()
		if y[0] == '#CHROM':
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
		chr = chr[3:]
	for i in range( 9, len(y) ):
		if i in select_id:
			tttt.append( y[i][0] )
			tttt.append( y[i][2] )
	info[ chr + ':' + pos ] = ''.join( tttt )
f.close()


def deal_title( t ):
	fr = {}
	chr, pos, ids = t[0], t[1], t[2:]
	for i in range( 0, len( ids ) ):
		tru_id, tru_hap = ids[i].split('_')[0], ids[i].split('_')[1]	
		if tru_id in d:
			fr[i+2] = [ tru_id, tru_hap ]
		else:
			print( 'Not.include' )
	return fr


def deal_one_row( y ):
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


fout = gzip.open( outfile, 'wt' )
f = gzip.open( infile, 'rt' )
title = f.readline().split()
re = deal_title( title ) # re is {}
#write title
ti = [ 'Chr', 'Pos' ]
for r in re:
	ti.append( '_'.join(re[r]) )
fout.write( '\t'.join( ti ) + '\n' )
#write rest
for line in f:
	y = line.split()
	chr, pos, results = y[0], y[1], y[2:]
	fw = deal_one_row( y )
	fout.write( '\t'.join( [ chr, pos ] + fw ) + '\n' )
f.close()
