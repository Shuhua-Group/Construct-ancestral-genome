#Each hap was paintered 10 times.
import sys,os,gzip

infile  = sys.argv[1] #chromopainter
inlist  = sys.argv[2] #id pop rank
outfile = sys.argv[3]
cutoff  = int(sys.argv[4]) #>=cutoff

d = {}
f = open( inlist )
for line in f:
	y = line.split()
	id, pop, rank = y[0], y[1], y[2]
	d[ rank ] = [ id, pop ]
f.close()


def deal_one_row( se, y, inf ):
	for k in range( 1, len(y) ):
		pop = d[ y[k] ][1]
		if k not in se:
			se[ k ] = { pop: 1 }
		else:
			if pop not in se[k]:
				se[ k ][ pop ] = 1
			else:
				se[ k ][ pop ] += 1
	return se

def deal_one_hap( u ):
	#u[0] was sample info row
	inf = u[0].split()
	se = {}
	for j in range( 1, len(u) ):
		y = u[j].split()
		se = deal_one_row( se, y, inf )

	for mem in se:
		p_sort = sorted( se[mem].items(), key = lambda x:x[1], reverse = True )
		if p_sort[0][1] >= cutoff:
			inf.append( p_sort[0][0] )
		else:
			inf.append( '.' )
	return inf
	

fout = gzip.open( outfile, 'wt' )
l = []
f = gzip.open( infile, 'rt' )
f.readline()
lines = f.readlines()
f.close()
n = len(lines)
for i in range( 0, int(n/11) ):
	start, end = 11 * i, 11 * i + 10
	fw = deal_one_hap( lines[ start: end + 1 ] )
	fout.write( '\t'.join( fw ) + '\n' )
fout.close()
