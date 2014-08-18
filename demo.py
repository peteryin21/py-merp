#!/bin/python

rs_dict = {}
with open('allrs',"r") as f:
	while True:
		l = f.readline()
		if l == '\n':
			break
		s = l.rstrip('\n')
		#print s
		if s not in rs_dict.keys():
			rs_dict[s] = []

	rs_set = set(rs_dict.keys())
	print len(rs_set)


written = {}
with open('pval' +"demo","w") as w:
	with open('bin/pval_file.txt',"r") as t:
		header = t.readline()
		w.write(header)
		while True:
			line = t.readline()
			if line == '\n' or line == '':
				break
			entry = line.rstrip('\n').split(' ')
			rs = entry[0]
			if rs in rs_set:
				written[rs] = []
				w.write(line)
		for rs in rs_set:
			if rs not in written.keys():
				print rs



