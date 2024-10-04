import sys,os,gzip

infile  = sys.argv[1]
poslist = sys.argv[2]
outfile = sys.argv[3]

#read 1
d = {}
k = 1
f = open( poslist )
for line in f:
	y = line.split()
	chr, pos = y[0], y[2]
	d[ k ] = [ chr, pos ]
	k += 1
f.close()

#read 2
l = []
f = gzip.open( infile, 'rt' )
for line in f:
	y = line.split()
	l.append( y )
f.close()

#write title
fout = gzip.open( outfile, 'wt' )
fout.write( 'Chr\tPos' )
for i in range( len(l) ):
	fout.write( '\t' + l[i][2] + '_' + l[i][1] )
fout.write( '\n' )

#write other
for i in range(1,len(d)+1):
	t = d[i]
	for j in range( len(l) ):
		t.append( l[j][2+i] )
	fout.write( '\t'.join( t ) + '\n' )
fout.close()
