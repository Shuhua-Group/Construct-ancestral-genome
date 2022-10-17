import sys,os,gzip,random

outfile        = sys.argv[1]
cutoff_freq    = float(sys.argv[2])
cutoff_length  = float(sys.argv[3])
infiles        = sys.argv[4:]

info = []
for infile in infiles:
	k = 0
	f = gzip.open( infile, 'rt' )
	f.readline()
	for line in f:
		y = line.split()
		chr, pos, genos = y[0], y[1], y[2:]
		if k + 1 > len(info):
			info.append( [ chr, pos, '|'.join(genos) ] )
			k += 1
		else:
			info[k][2] = info[k][2] + '|' + '|'.join(genos)
			k += 1
	f.close()


def deal_one_row( point, n, slgeno ):
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


def judge( li ):
	kl = 0
	for i in range( len(li) ):
		if li[i] != '.':
			kl += 1
	#print( '%.6f'%(float(kl) / len(li) ) )
	if cutoff_freq <= float(kl) / len(li):
		return 'ok'
	else:
		return 'no'

fout = gzip.open( outfile, 'wt' )
n = len(info[0][2].split('|'))
point = random.randint( 0, n-1 )
start = '.'
for i in range( len(info) ):
	slgeno = info[i][2].split('|')
	res = judge( slgeno )
	if res != 'ok':
		point = random.randint( 0, n-1 )
		fout.write( info[i][0] + '\t' + info[i][1] + '\t' + '.' + '\n' )
		start = '.'
		continue
	if slgeno[point] == '.':
		fw, point = deal_one_row( point, n, slgeno )
		if fw == '.':
			fout.write( info[i][0] + '\t' + info[i][1] + '\t' + '.' + '\n' )
			start = '.'
			#start, ind = '.', '.'
		else:
			fout.write( info[i][0] + '\t' + info[i][1] + '\t' + fw  + '\n' )
			start = info[i][1]
			#start, ind = info[i][1], fw[2:]
	else:
		result = ''
		if start == '.':
			result = 'ok'
		else:
			if int(info[i][1]) + 1 - int(start) > cutoff_length:
				result = 'no'
		if result == 'no':
			fw, point = deal_one_row( point, n, slgeno )
			if fw == '.':
				fout.write( info[i][0] + '\t' + info[i][1] + '\t' + fw  + '\n' )
				start = '.'
				#start, ind = '.', '.'
			else:
				fout.write( info[i][0] + '\t' + info[i][1] + '\t' + fw  + '\n' )
				start = info[i][1]
				#start, ind = info[i][1], fw[2:]
		else:
			fout.write( info[i][0] + '\t' + info[i][1] + '\t' + slgeno[point] + '\n' )
fout.close()
		
