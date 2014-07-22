###rs9282541 not in 1000genomes, what to do when this happens? likely a common problem

####ADD HEADER COLumn FOR NEW ALLELE
#!/bin/python
import sys
import pdb
import os
trait_list = []

if len(sys.argv) == 1:
	print "Usage: python update.py [trait] [trait2]..."
	sys.exit()
for trait in sys.argv[1:]:
	trait_list.append(trait)
	

	
allele_file = "data/1000_genomes"

negative_list = ["decrease","lower","shorter"]
nucleotides = ['A','T','G','C']
allele_handle = file(allele_file, "r")
allele_lines = allele_handle.readlines()
allele_dict = {}

def main(trait):
	def get_complement(nuc):
		if nuc == "A":
			return "T"
		if nuc == "T":
			return "A"
		if nuc == "G":
			return "C"
		if nuc == "C":
			return "G"


	handle = file(trait,"r")
	lines = handle.readlines()
	snp_dict = {}

	write_handle = file(trait+"_update", "w")
	#####MANUALLY CHANGE HEADER TO INCLUDE NON RISK ALLELE COLUMN
	header = lines[0]
	write_handle.write(header)
	for line in lines[1:]:
		entry = line.rstrip('\n').split('\t')
		risk_allele = entry[6]
		snp = entry[5]
		unit = entry[4]
		snp_dict[snp] = risk_allele

	handle.close()
	non_risk_dict = {}
	for line in allele_lines:
		entry = line.rstrip('\n').split('\t')
		snp = entry[5]

		if snp in snp_dict.keys():
			allele_one = entry[2]
			allele_two = entry[3]
			risk = snp_dict[snp]
			if snp == "rs560887":
				pdb.set_trace()
			if risk == allele_one:
				non_risk_dict[snp] = allele_two
			elif risk == allele_two:
				non_risk_dict[snp] = allele_one
			else:
				###if neither allele is the risk, take complements and work from there
				complement = get_complement(risk)
				if complement == allele_one:
					non_risk_dict[snp] = get_complement(allele_two)
				if complement == allele_two:
					non_risk_dict[snp] = get_complement(allele_one)

	handle = file(trait,"r")
	lines = handle.readlines()
	not_in_genomes = ""
	for line in lines[1:]:
		entry = line.rstrip('\n').split('\t')
		risk_allele = entry[6]
		snp = entry[5]
		unit = entry[4]
		snp_dict[snp] = risk_allele
		try:
			allele_entry = non_risk_dict[snp]
		except:
			allele_entry = "CHECK_dbSNP_for_" + snp
			not_in_genomes = not_in_genomes + snp + '\n'

		written = False
		for word in negative_list:
			if word in unit:
				entry_unit = -(float(entry[3]))
				entry_new = str(entry_unit)
				write_handle.write(entry[0] + '\t' + entry[1] + '\t' + entry[2] + '\t' + entry_new + '\t' + entry[4]
					+ '\t' + entry[5]+ '\t' + allele_entry + '\t' + entry[6]+ '\t' + entry[7]+ '\t' + entry[8] + '\n' )
				written = True
		if written == False:
			write_handle.write(entry[0] + '\t' + entry[1] + '\t' + entry[2] + '\t' + entry[3] + '\t' + entry[4]	+ '\t' + entry[5]+ '\t' + allele_entry + '\t' + entry[6]+ '\t' + entry[7]+ '\t' + entry[8] + '\n' )
	if not_in_genomes != "":
		print not_in_genomes + "Above are SNPs that are not in 1000 genomes database for " + trait + ". Manually enter the non-risk allele for these SNPs"

for trait in trait_list:
	main(trait)

##loop back through trait file and add beta sign changing along with rewriting all lines 
##combine gwas pull?
