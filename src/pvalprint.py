#!/bin/python

rs = []
with open('bmd_list.txt',"r") as f:
	while True:
		line = f.readline()
		if line == "":
			break
		entry = line.rstrip('\n')
		rs.append(entry)

	rs_set = set(rs)

with open('pval_file.txt','r') as p:
	while True:
		line = p.readline()
		if line == "":
			break
		entry = line.rstrip('\n').split(' ')
		if entry[0] in rs_set:
			print line.strip()



