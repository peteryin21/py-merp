##add filter algorithm componenets here, combine
#!/bin/python

#Paramters
#Filtering for associations in allmetabolic pval file
#pmax 1:if weaker association, if a SNP has more than threshold1 with p<pmax1, SNP excluded.
pmax1 = 0.01
threshold1 = 3 #4 or more are taken out
#pmax 2: if there is even one strong association of pmax2, SNP excluded.
pmax2 = 0.001
threshold2 = 0 # None can be less than .001

###LD filtering####
#SNPs with R^2 value > rsq_threshold are clustered together, with only lead SNP passing on.
rsq_threshold = 0.05


##Pval iteration###
#If number of total pmax1 violations exceeds max_fraction * total number of tests, then cut SNPS until below max_fraction of tests.
max_fraction = 0.05
#0.05



import sys
import pdb
import os
import requests
include_list = []
exclude_list = []

num_args = len(sys.argv) - 1
try:
	trait_file = sys.argv[1]
	include_file = sys.argv[2]
	exclude_file = sys.argv[3]
	
except:
	print "Usage: python filter.py [trait_file] [include_traits_file] [exclude_file]"
	sys.exit()
#Read include list file and exclude list file(pval headers)	

###GETTING LD FILE FROM LD SNAP BROAD PROXY USING REQUESTS######
rs_list = []
with open(trait_file,"r") as trait:
	lines = trait.readlines()
	for line in lines[1:]:
		entry = line.rstrip('\n').split('\t')
		rs_list.append(entry[5])
rs = ('\n\t').join(rs_list)
#rs = '"""'+'\n' +rs+ '\n' + '"""'
rs =  '\n\t' +rs 
#pdb.set_trace()
print rs

# rs = """
# 	rs9627183
# 	rs11134178
# 	rs12915721
# 	rs2157697
# 	rs4011946
# 	rs6501530
# 	rs12301774
# 	rs2594278
# 	rs3792452 
# 	rs1065758
# 	rs12360508
# 	rs16847570
# 	rs175126
# 	rs933771
# 	rs1381795
# 	rs3845659
# 	rs2382075
# 	rs16959263
# 	rs1005324
# 	rs2039430
# 	rs4853259
# 	rs12651081
# 	rs11999224
# 	rs17755054
# 	rs11652864
# 	rs2653306
# 	rs7791083
# 	rs10217716
# 	rs6580967"""
data_1000 = {
	"searchPairwise": "true",
	"snpList": rs,
	"hapMapRelease": "onekgpilot",
	"hapMapPanel": "CEU",
	"RSquaredLimit": "0",
	"distanceLimit": "500000",
	"downloadType": "Browser",
	"arrayFilter": "query",
	"columnList[]": "DP",
	"columnList[]": "GP",
	"columnList[]": "AM",
	"submit":"search"
}
##backup ld search in case not in 1000genomes build
data_hm22 = {
	"searchPairwise": "true",
	"snpList": rs,
	"hapMapRelease": "rel22",
	"hapMapPanel": "CEU",
	"RSquaredLimit": "0",
	"distanceLimit": "500000",
	"downloadType": "Browser",
	"arrayFilter": "query",
	"columnList[]": "DP",
	"columnList[]": "GP",
	"columnList[]": "AM",
	"submit":"search"
}
r_1000 = requests.post('http://www.broadinstitute.org/mpg/snap/ldsearch.php', data=data_1000)
ld_result = r_1000.text.split('\n')

r_hm22 = requests.post('http://www.broadinstitute.org/mpg/snap/ldsearch.php', data=data_hm22)
ld_result_hm22= r_hm22.text.split('\n')
##create list of lines of ld_result

# pdb.set_trace()

#LD_file = "bin/1LD.txt"

nhgri_file = "data/finalgwas.txt"
pval_file = "data/v3abr_allmetabolic_pvals_v2.txt"
replacement_pval_snps = "data/replacement_pval_snps.txt"


with open(include_file,"r") as include:
	lines = include.readlines()
	for line in lines:
		entry = line.rstrip('\n')
		include_list.append(entry)

with open(exclude_file,"r") as excluded_headers:
	lines = excluded_headers.readlines()
	for line in lines:
		entry = line.rstrip('\n')
		exclude_list.append(entry)



#Add all diseases we want to exlcude from NHGRI filtering to include_list
disease_list_handle = file("./data/disease_list.txt","r")
disease_line = disease_list_handle.readlines()
for line in disease_line:
	entry = line.rstrip('\n').split('\t')
	disease = entry[0]
	include_list.append(disease)


####NHGRI PORTION#######
##Creates dictionary nhgri_dict mapping SNP rs# to list of traits it's associated with##
##Creates dictionary snp_nhgri_dict mapping snp to true if it passes nhgri test and false if not##

trait_handle = file(trait_file,"r")
trait_lines = trait_handle.readlines()
snp_list = []
for line in trait_lines[1:]:
	line_list = line.rstrip('\n').split('\t')
	snp = line_list[5]
	if snp not in snp_list:
		snp_list.append(snp)
trait_handle.close()
print str(len(snp_list)) + " is number of snps"

nhgri_handle = file(nhgri_file, "r")
nhgri_lines = nhgri_handle.readlines()
header = trait_lines[0]

nhgri_dict = {}
for line in nhgri_lines[1:]:
	if line == '\n':
		break
	mod_line = line.rstrip('\n').split('\t')
	rs = mod_line[21]
	trait = mod_line[7]
	if rs in snp_list:
		if rs in nhgri_dict.keys():
			nhgri_dict[rs].append(trait)
		else:
			nhgri_dict[rs] = [trait]
nhgri_handle.close()

def trait_included(trait,include_list):
	for e in include_list:
		if e in trait:
			return True
	return False



def nhgri_test(snp,include_list):
	for trait in nhgri_dict[snp]:
		if trait_included(trait,include_list):
			pass
		else:
			print snp + " is also associated with " + trait

			return False

	return True
					

		


snp_nhgri_dict = {}
for snp in snp_list:
	snp_nhgri_dict[snp] = nhgri_test(snp,include_list)	

####NHGRI PORTION END#######
		
######PVAL PORTION START########
pval_handle = file(pval_file, "r")





dict_snp = {}
pval_lines = pval_handle.readlines()
header = pval_lines[0]
header_split = header.rstrip('\n').split(' ')
columns = len(header_split) - 3
correction_num = 0
for p in header_split:
	for excluded_trait in exclude_list:
		if excluded_trait in p:
			print str(columns) + "before"
			columns = columns - 1
			print str(columns) + "now"
			correction_num += 1






# replace_handle = file(replacement_pval_snps, 'r')
# replace_lines = replace_handle.readlines()
# replace_dict = {}
# for line in replace_lines[1:]:
# 	entry = line.rstrip('\n').split('\t')
# 	orig = entry[0]
# 	replacement = entry[1]
# 	if orig in no_pval_snps.keys():
# 		if orig not in replace_dict.keys():
# 			replace_dict[orig] = replacement
# 		else:
# 			print orig + " is repeated in replacement text file"
# 			pass



for line in pval_lines[1:]:
	mod_line = line.rstrip('\n').split(' ')
	rs = mod_line[0]
	count_p1 = 0
	count_p2 = 0
	# if rs == "rs2366858":
	# 	pdb.set_trace()
#################################
	#make height all = 1 to remove from 
	mod_line[9] = 1
###################################
	#Will be more efficient to move this in next if statement? For v2
	for e in mod_line[3:]:

		if float(e) <= pmax1:
			count_p1 = count_p1 + 1
		if float(e) <= pmax2:
			count_p2 = count_p2 + 1
	count_p1 = count_p1 - correction_num
	count_p2 = count_p2	- correction_num

	if rs not in dict_snp.keys() and rs in snp_list: #only add to dict if rs in trait file
		#first compare true count to modified threshold then change count to modified count if trait_related is true
		if count_p1 <= threshold1:
			if count_p2 <= threshold2:
				#if trait in pval file, subtract 1 from count of number of sig associations
				dict_snp[rs] = [True, count_p1] # Record number of associations wiht p<0.05 afte rfiltering for less than 4
			else:
				dict_snp[rs] = [False, count_p1]
		else:
			dict_snp[rs] = [False, count_p1]
	

	# elif rs in replace_dict.keys():
	# 	if count_p1 <= threshold1:
	# 		if count_p2 <= threshold2:
	# 			#if trait in pval file, subtract 1 from count of number of sig associations
	# 			dict_snp[rs] = [True, count_p1] # Record number of associations wiht p<0.05 afte rfiltering for less than 4
	# 		else:
	# 			dict_snp[rs] = [False, count_p1]
	# 	else:
	# 		dict_snp[rs] = [False, count_p1]
pval_handle.close()

###IF NOT IN PVAL### remove from snp_list for now, later change so replaces or something else!!
no_pval_snps = {}
for snp in snp_list:
	if snp not in dict_snp.keys():
		no_pval_snps[snp] = []
		# snp_list.remove(snp)
###calculate total number of significant assoc (<0.05)

num_sig = 0
trait_handle = file(trait_file,"r")
trait_lines = trait_handle.readlines()
#pdb.set_trace()
for line in trait_lines[1:]:
	line_list = line.rstrip('\n').split('\t')
	if len(line_list) > 5:
		snp = line_list[5]
	else:
		pass
	#pdb.set_trace() 
	if snp in dict_snp:
		#Create new dict with only trait files
		#abridged_dict[snp] = dict_snp[snp] # should be unnecessary

		if dict_snp[snp][0] == True: 
			num_sig = num_sig + dict_snp[snp][1]
		#if trait in pval file, subtract 1 from num_sig
			#if trait_related == True:
				#num_sig = num_sig - 1
trait_handle.close()

####PVAL PORTION END #####
####LD PORTION BEGIN######

trait_handle = file(trait_file,"r")
pval_dict = {}
repeated_snps = []
trait_lines = trait_handle.readlines()
for line in trait_lines[1:]:
	line_list = line.rstrip('\n').split('\t')
	snp = line_list[5]
	p = line_list[2]
	if snp not in pval_dict:
		pval_dict[snp] = p
	else:
		p_new = p
		#replace p only if more significant
		if float(p_new) < float(pval_dict[snp]):
			pval_dict[snp] = p_new
		if snp not in repeated_snps:
			repeated_snps.append(snp)
trait_handle.close()

#Go through and rewrite traitfile without repeating snps
trait_handle = file(trait_file,"r")
trait_lines = trait_handle.readlines()
header = trait_lines[0]
abridged_trait_handle = file("./"+trait_file + "_abr_temp","w")
abridged_trait_handle.write(header)
for line in trait_lines[1:]:
	line_list = line.split('\t')
	snp = line_list[5]
	p = line_list[2]
	pid = line_list[9]
	if snp not in repeated_snps:
		abridged_trait_handle.write(line)
	else:
		if float(p) == float(pval_dict[snp]):
			abridged_trait_handle.write(line)

abridged_trait_handle.close()
trait_handle.close()

cluster_dict = {}
index_dict = {}

#ld_handle = file(LD_file,"r")
#ld_lines = ld_handle.readlines()


# def ld_cluster(snp,proxy,line_list,cluster_index,not_in_ld,ld_assoc_dict,non_assoc_snps):
# if "WARNING" in proxy: #and "query" in proxy #and "not" in proxy and "in" in proxy:
# 	if "Query snp not in" in line_list[2]: 
# 		if snp not in not_in_ld:
# 			#
# 			not_in_ld.append(snp)
# 		pass # should go back to line in ld_lines
# 	elif "No matching proxy snps found" in line_list[2]:
# 		#This means that snp is it's own cluster, treat as cluster with only one element
# 		if snp not in not_in_ld and snp not in non_assoc_snps:
# 			cluster_dict[cluster_index] = [snp]
# 			index_dict[snp] = cluster_index
# 			cluster_index +=1
# 			non_assoc_snps.append(snp)

# elif "rs" in proxy:
# 	rsq = line_list[3]

# 	#Create dict of association of tuple to rsq
# 	ld_assoc_dict[(snp,proxy)] = rsq

# 	if float(rsq) >= rsq_threshold:
# 		if snp not in index_dict.keys() and proxy not in index_dict.keys():
# 			cluster_dict[cluster_index] = [snp,proxy]
# 			index_dict[snp] = cluster_index
# 			index_dict[proxy] = cluster_index
# 			cluster_index += 1
# 		elif snp in index_dict.keys() and proxy not in index_dict.keys():
# 			snp_index = index_dict[snp]
# 			cluster_dict[snp_index].append(proxy)
# 			index_dict[proxy] = snp_index

# 		elif snp not in index_dict.keys() and proxy in index_dict.keys():
# 			proxy_index = index_dict[proxy]
# 			cluster_dict[proxy_index].append(snp)
# 			index_dict[snp] = proxy_index
# 		elif snp in index_dict.keys() and proxy in index_dict.keys():
# 			snp_index = index_dict[snp]
# 			proxy_index = index_dict[proxy]
# 			if snp_index != proxy_index:
# 				#change b index to as
# 				index_dict[proxy] = snp_index
# 				for e in cluster_dict[proxy_index]:
# 					index_dict[e] = snp_index
# 				#combine lists 
# 				cluster_dict[snp_index] = cluster_dict[snp_index] + cluster_dict[proxy_index]
# 				#delete entry for b
# 				del cluster_dict[proxy_index]
# 			else:
# 				pass

# pdb.set_trace()

cluster_index = 1
not_in_ld = []
ld_assoc_dict = {}
non_assoc_snps = []
for line in ld_result[1:]:
	line_list = line.split('\t')
	if len(line_list) != 1:
		snp = line_list[0]

		proxy = line_list[1]
		if "WARNING" in proxy: #and "query" in proxy #and "not" in proxy and "in" in proxy:
			if "Query snp not in" in line_list[2]: 
				if snp not in not_in_ld:
			#
					not_in_ld.append(snp)
				pass # should go back to line in ld_lines
			elif "No matching proxy snps found" in line_list[2]:
		#This means that snp is it's own cluster, treat as cluster with only one element
				if snp not in not_in_ld and snp not in non_assoc_snps:
					cluster_dict[cluster_index] = [snp]
					index_dict[snp] = cluster_index
					cluster_index +=1
					non_assoc_snps.append(snp)

		elif "rs" in proxy:
			rsq = line_list[3]

			#Create dict of association of tuple to rsq
			ld_assoc_dict[(snp,proxy)] = rsq

			if float(rsq) >= rsq_threshold:
				if snp not in index_dict.keys() and proxy not in index_dict.keys():
					cluster_dict[cluster_index] = [snp,proxy]
					index_dict[snp] = cluster_index
					index_dict[proxy] = cluster_index
					cluster_index += 1
				elif snp in index_dict.keys() and proxy not in index_dict.keys():
					snp_index = index_dict[snp]
					cluster_dict[snp_index].append(proxy)
					index_dict[proxy] = snp_index

				elif snp not in index_dict.keys() and proxy in index_dict.keys():
					proxy_index = index_dict[proxy]
					cluster_dict[proxy_index].append(snp)
					index_dict[snp] = proxy_index
				elif snp in index_dict.keys() and proxy in index_dict.keys():
					snp_index = index_dict[snp]
					proxy_index = index_dict[proxy]
					if snp_index != proxy_index:
						#change b index to as
						index_dict[proxy] = snp_index
						for e in cluster_dict[proxy_index]:
							index_dict[e] = snp_index
						#combine lists 
						cluster_dict[snp_index] = cluster_dict[snp_index] + cluster_dict[proxy_index]
						#delete entry for b
						del cluster_dict[proxy_index]
					else:
						pass
		




###HERE GO THROUGH SAME THING FOR LD-result hm2 except check for only not inld

				
if len(not_in_ld)>=1:
	for line in ld_result_hm22[1:]:
		line_list = line.split('\t')
		if len(line_list) != 1:
			snp = line_list[0]
			# if snp == "rs2366858":
			# 	pdb.set_trace()
			proxy = line_list[1]
			if snp in not_in_ld:
				if "WARNING" in proxy: #and "query" in proxy #and "not" in proxy and "in" in proxy:
					if "Query snp not in" in line_list[2]: 
						if snp not in not_in_ld:
					#
							not_in_ld.append(snp)
						pass # should go back to line in ld_lines
					elif "No matching proxy snps found" in line_list[2]:
				#This means that snp is it's own cluster, treat as cluster with only one element
						if snp not in not_in_ld and snp not in non_assoc_snps:
							cluster_dict[cluster_index] = [snp]
							index_dict[snp] = cluster_index
							cluster_index +=1
							non_assoc_snps.append(snp)

				elif "rs" in proxy:
					rsq = line_list[3]

					#Create dict of association of tuple to rsq
					ld_assoc_dict[(snp,proxy)] = rsq

					if float(rsq) >= rsq_threshold:
						if snp not in index_dict.keys() and proxy not in index_dict.keys():
							cluster_dict[cluster_index] = [snp,proxy]
							index_dict[snp] = cluster_index
							index_dict[proxy] = cluster_index
							cluster_index += 1
						elif snp in index_dict.keys() and proxy not in index_dict.keys():
							snp_index = index_dict[snp]
							cluster_dict[snp_index].append(proxy)
							index_dict[proxy] = snp_index

						elif snp not in index_dict.keys() and proxy in index_dict.keys():
							proxy_index = index_dict[proxy]
							cluster_dict[proxy_index].append(snp)
							index_dict[snp] = proxy_index
						elif snp in index_dict.keys() and proxy in index_dict.keys():
							snp_index = index_dict[snp]
							proxy_index = index_dict[proxy]
							if snp_index != proxy_index:
								#change b index to as
								index_dict[proxy] = snp_index
								for e in cluster_dict[proxy_index]:
									index_dict[e] = snp_index
								#combine lists 
								cluster_dict[snp_index] = cluster_dict[snp_index] + cluster_dict[proxy_index]
								#delete entry for b
								del cluster_dict[proxy_index]
							else:
								pass
#pdb.set_trace()		
#now take out snps from notinld that were found in hm22
#
for snp in not_in_ld:
	# if snp == "rs2366858":
	# 		pdb.set_trace()
	for key in index_dict.keys():
		if snp==key: #test
			not_in_ld.remove(snp)
# for line in ld_lines[1:]:
# 	line_list = line.rstrip('\n').split('\t')
# 	snp = line_list[0]
# 	proxy = line_list[1]
# 	if "WARNING" in proxy: #and "query" in proxy #and "not" in proxy and "in" in proxy:
# 		if "Query snp not in" in line_list[2]: 
# 			if snp not in not_in_ld:
# 				not_in_ld.append(snp)
# 			pass # should go back to line in ld_lines
# 		elif "No matching proxy snps found" in line_list[2]:
# 			#This means that snp is it's own cluster, treat as cluster with only one element
# 			if snp not in not_in_ld and snp not in non_assoc_snps:
# 				cluster_dict[cluster_index] = [snp]
# 				index_dict[snp] = cluster_index
# 				cluster_index +=1
# 				non_assoc_snps.append(snp)

# 	elif "rs" in proxy:
# 		rsq = line_list[3]

# 		#Create dict of association of tuple to rsq
# 		ld_assoc_dict[(snp,proxy)] = rsq

# 		if float(rsq) >= rsq_threshold:
# 			if snp not in index_dict.keys() and proxy not in index_dict.keys():
# 				cluster_dict[cluster_index] = [snp,proxy]
# 				index_dict[snp] = cluster_index
# 				index_dict[proxy] = cluster_index
# 				cluster_index += 1
# 			elif snp in index_dict.keys() and proxy not in index_dict.keys():
# 				snp_index = index_dict[snp]
# 				cluster_dict[snp_index].append(proxy)
# 				index_dict[proxy] = snp_index

# 			elif snp not in index_dict.keys() and proxy in index_dict.keys():
# 				proxy_index = index_dict[proxy]
# 				cluster_dict[proxy_index].append(snp)
# 				index_dict[snp] = proxy_index
# 			elif snp in index_dict.keys() and proxy in index_dict.keys():
# 				snp_index = index_dict[snp]
# 				proxy_index = index_dict[proxy]
# 				if snp_index != proxy_index:
# 					#change b index to as
# 					index_dict[proxy] = snp_index
# 					for e in cluster_dict[proxy_index]:
# 						index_dict[e] = snp_index
# 					#combine lists 
# 					cluster_dict[snp_index] = cluster_dict[snp_index] + cluster_dict[proxy_index]
# 					#delete entry for b
# 					del cluster_dict[proxy_index]
# 				else:
# 					pass
print "check cluster_dict"
pdb.set_trace()
def most_sig(snp_list):
	sig_snp = snp_list[0]
	for e in snp_list:
		if float(pval_dict[e]) <= float(pval_dict[sig_snp]):
			sig_snp = e
	return sig_snp

# rsq_nhgri_threshold = 0.8
# rsq_nhgri_lower_threshold = 0.1
#pdb.set_trace()
abr_trait_handle = file(trait_file+"_abr_temp","r")
trait_lines = abr_trait_handle.readlines()
header = trait_lines[0]
cluster_status = {}
most_sig_cluster = {}

for line in trait_lines[1:]:
	line_list = line.split('\t')
	snp = line_list[5]
	# if snp == "rs2366858":
	# 	pdb.set_trace()
	p = line_list[2]
	if snp in not_in_ld:
		pass
	elif snp in index_dict.keys():
		index = index_dict[snp]
		assoc_list = cluster_dict.get(index)
		sig_rs = most_sig(assoc_list)

		most_sig_cluster[index] = sig_rs
		#Case where most sig snp not found in pval file
		# case = False
		# pdb.set_trace()
		# while case == False:
		# 	if sig_rs not in dict_snp.keys():
		# 		sig_rs = most_sig(assoc_list.remove(sig_rs))
		# 	if sig_rs in dict_snp.keys():
		# 		case == True

		# if case == False:
		# 	pass






		print sig_rs + " is the most sig out of cluster " + str(index)
abr_trait_handle.close()
pdb.set_trace()
for key in most_sig_cluster.keys():
	snp = most_sig_cluster[key]
	# if snp == "rs2366858":
	# 	pdb.set_trace()
	index = index_dict[snp]
	assoc_list = cluster_dict.get(index)
	assoc_list.remove(snp)
	#nhgri test
	if snp_nhgri_dict[snp] == True:
		##pval test\
		if snp in dict_snp.keys():
			if dict_snp[snp][0] == True: 
				
				cluster_status[snp] = True

				###1/2/14. I think this section is the part that makes sure that if a SNP is LD assoc
				#with a SNP that fails nghri filter, that it is kicked out, though we use a different
				#threshold for this, higher one. Should add this back in, v2, also check for pval violations with assoc snps
				#
				###WE ONLY CARE ABOUT LEAD SNP IN EACH ONE, EVEN IF PASSES, SHOULD WE INLCUDE THINGS IN LD?
				#has passed both nhgri and pval test
				#now check to see if rsq between lead snp and each rs in cluster is greather than nhgri threshold of 0.1, 
				#more strict? Commenting out for now till understand i think its not necessary
				# for rs in assoc_list:
				# 	rsq = ld_assoc_dict[(snp,rs)]

				# 	if rsq >= rsq_nhgri_threshold:
				# 		if snp_nhgri_dict[rs] == False:
				# 			cluster_status[snp] = False
				# 			break
				# 		else: 
				# 			pass
				# 			#elif rsq >= rsq_nhgri_lower_threshold:
				# 			#	if snp_nhgri_dict[rs] == False:
				# 			#		cluster_status[snp] = "IDK"
				# 			#		print "ambigious for" + snp
				# 			#		pass
				# 			#	else: 
				# 			cluster_status[snp] = True
				# 	else:
				# 		cluster_status[snp] = True
			else:
				cluster_status[snp] = False
		else:
			cluster_status[snp] = False
	else:
		print snp + 'not in pval file, fix this'
		cluster_status[snp] = False

#######pval violation >0.05 by chance algorithm#####
def dict_count(d):
	count = 0
	for key in d:
		#If passes conditions count it
		if d[key][0] == True:
			count += 1
	print count
	return count
#Help function takes dictionary and cutoff and turns all keys with num_sig greater/= to cutoff to false
#Then returns number of sig
def dict_purge_count(d,cluster_d,cut,old_num_sig):
	#pdb.set_trace()
	return_sig = 0
	for key in d:
		if key in cluster_d:
			if d[key][1] >= cut:
				cluster_d[key] = False 
 
	#get new number of sig
##currently ignore if not in pval file	
	for snp in dict_snp.keys():
		
		if snp in d: #unsure if dict_snp or dict
			if snp in cluster_d:
				if cluster_d[snp] == True: 
					return_sig = return_sig + d[snp][1]

	# for snp in snp_list:
	# 	if snp == "rs4826508":
	# 		pdb.set_trace()
	# 	if key in d: #unsure if dict_snp or dict
	# 		if key in cluster_d:
	# 			if cluster_d[key] == True: 
	# 				return_sig = return_sig + d[snp][1]

	print str(return_sig) + "is new sig after round changing cut"
	return return_sig


threshold_met = False
new_num_sig = num_sig
cut = threshold1


while threshold_met == False:
	#num_snps = dict_count(abridged_dict)

	####JUST LOOK AT ELIGIBLE SNPS
	num_snps = 0
	for key in cluster_status.keys():
		if cluster_status[key] == True:
			num_snps +=1

	num_tests = num_snps * columns
	#new_num_sig = num_sig
	if float(new_num_sig) <= max_fraction * num_tests:
		threshold_met = True
	#pdb.set_trace()
	if float(new_num_sig) > max_fraction * num_tests:
		#ok??
		#pdb.set_trace()
		old_num_sig = new_num_sig
		new_num_sig = dict_purge_count(dict_snp,cluster_status,cut,old_num_sig)
		cut = cut - 1
		if cut < 0:
			print "Fatal error: cut below 0, consider having a higher cutoff threshold of violations or higher max_fraction "
			sys.exit()
#############################
# print "see if rs2366858 is in dictsnp"
# pdb.set_trace()

snps_to_write = []
non_assoc_write_temp = []
for key in cluster_status.keys():
	if cluster_status[key] == True:
		snps_to_write.append(key)
# path="./traitFiles/final_trait_files"
# if not os.path.exists(path):
#     os.makedirs(path)
###cont this later
abr_trait_handle = file(trait_file+"_abr_temp","r")
trait_lines = abr_trait_handle.readlines()
header = trait_lines[0]
updated_handle = file(trait_file +"filtered2", "w")
updated_handle.write(header)
for line in trait_lines[1:]:
	line_list = line.split('\t')
	snp = line_list[5]
	if snp in snps_to_write:
		updated_handle.write(line)
	# if snp in non_assoc_snps:
	# 	updated_handle.write("**Need to check if violating, repeated in file??**" + line)

#DGAF about them snps not in ld anymore
# for snp in not_in_ld:
# 	updated_handle.write(snp + " ---Not in LD Data---" + '\n')

os.remove(trait_file+"_abr_temp")


			#if snp in pval_dict:
		#	else: 
		#		print snp + " not counted because of nhgri filter"
		##		pass
		#3else:
		#	print snp + " does not make the LD cut"
	#else:
	#	final_trait_handle.write(line)
	
#add snps not in LD database at end
#for e in not_in_ld:
	#print "YO"
	#final_trait_handle.write("+" + e)

abr_trait_handle.close()

#print "E"
#dict 1 map 1,2,3 to a list of associated snps
#dict 2 map snps to their associated index, 1 2 3 etc.

#A B
#C D
#A E

#If A and B are both in dict 1 and significant assoc, change b index to A's in dict 2, add B and all clusters with B
# to index of A's cluster list and remove index of B's entry of cluster lists.
#If A or B is in dict 1 but the other is not, simply add the new snp to the cluster index in dict 1 and index in dict 2
#If neither A or B are in dict 1, create a new cluster and index in dict 1 and ref in dict 2
