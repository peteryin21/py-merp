#!/bin/python
import math
import pdb
import sys
import os

from scipy import stats

try:
    trait_file= sys.argv[1]
    disease_file = sys.argv[2]
except:
    print "Usage: python v3mr_calc_fin.py [path/to/traitfile] [path/to/diseasefile]"
    sys.exit()

######TRAIT FILE PARAMETERS######
beta_index = 3
rs_index = 5
allele1_index = 6
riskallele_index = 7


####DISEASE FILE PARAMTERS#######
dis_rs_index = 0
dis_allele1_index = 2
dis_riskallele_index = 4
dis_lnOR_index = 11
dis_lnse_index = 12



trait_handle = file(trait_file,"r")
dict_trait = {}
trait_content = trait_handle.readlines()
#pdb.set_trace()
#####


for entry in trait_content[1:]:
    #pdb.set_trace()
    if entry == '\n' or entry == ' ':
        break
    entry = entry.rstrip('\n').split('\t')
    beta = entry[beta_index]
    rsid = entry[rs_index]
    allele1 = entry[allele1_index]
    riskallele = entry[riskallele_index]
    #SNPrs# for now
#    pdb.set_trace()
    if not(rsid in dict_trait):
        dict_trait[rsid] = []
        dict_trait[rsid].append(beta)
        dict_trait[rsid].append(allele1)
        dict_trait[rsid].append(riskallele)
    else:
        #shouldn't happen
        pass
        #dict_trait[entry[0]].append( (entry[2]) )

 #Change to actual disease file input, 2nd input
disease_handle = file(disease_file,"r")
disease_content = disease_handle.readlines()
#pdb.set_trace()

#more efficient way
#Loop through disease file and if SNP ins trait file, add info

for entry in disease_content[1:]:
    entry = entry.rstrip('\n').split(' ')
 #   effallele = entry[4]
    if entry[dis_rs_index] in dict_trait.keys():
        #print entry[0]
        effallele = entry[dis_riskallele_index]

        dict_trait[entry[dis_rs_index]].append(entry[dis_lnOR_index])
        dict_trait[entry[dis_rs_index]].append(entry[dis_lnse_index])
        #eventually calculate from 95ci
    #now dict should be SNP: [beta_trait, allele1, riskallele, beta_mi,SE]
      # pdb.set_trace()
        if effallele == dict_trait[entry[dis_rs_index]][2]:
            pass#print "Aligned!"
        elif effallele == dict_trait[entry[dis_rs_index]][1]:
            #print "FLIP THE BETA"
            dict_trait[entry[dis_rs_index]][0]= -(float(dict_trait[entry[dis_rs_index]][0]))
            #If eff allele is not present in either allele 1 or risk, fix strand problem
        else:
            if effallele == 'A':
                if 'T' == dict_trait[entry[dis_rs_index]][2]:
                    pass#print "Aligned finally"
                elif 'T' == dict_trait[entry[dis_rs_index]][1]:
                    #print "Flip the BETA FINALLY"
                    dict_trait[entry[dis_rs_index]][0]= -(float(dict_trait[entry[dis_rs_index]][0]))
                else:
                    pass#print "wut" + effallele + '\t' + dict_trait[entry[0]][1] + dict_trait[entry[0]][2] 
            if effallele == 'T':
                if 'A' == dict_trait[entry[dis_rs_index]][2]:
                    pass#print "Aligned finally"
                elif 'A' == dict_trait[entry[dis_rs_index]][1]:
                    #print "Flip the BETA FINALLY"
                    dict_trait[entry[dis_rs_index]][0]= -(float(dict_trait[entry[dis_rs_index]][0]))
                else:
                    pass#print "wut" + effallele + '\t' + dict_trait[entry[0]][1] + dict_trait[entry[0]][2] 
            if effallele == 'C':
                if 'G' == dict_trait[entry[dis_rs_index]][2]:
                    pass#print "Aligned finally"
                elif 'G' == dict_trait[entry[dis_rs_index]][1]:
                    #print "Flip the BETA FINALLY"
                    dict_trait[entry[dis_rs_index]][0]= -(float(dict_trait[entry[dis_rs_index]][0]))
                else:
                    pass#print "wut" + effallele + '\t' + dict_trait[entry[0]][1] + dict_trait[entry[0]][2] 
            if effallele == 'G':
                if 'C' == dict_trait[entry[dis_rs_index]][2]:
                    pass#print "Aligned finally"
                elif 'C' == dict_trait[entry[dis_rs_index]][1]:
                    #print "Flip the BETA FINALLY"
                    dict_trait[entry[dis_rs_index]][0]= -(float(dict_trait[entry[dis_rs_index]][0]))
                else:
                    pass#print "wut" + effallele + '\t' + dict_trait[entry[0]][1] + dict_trait[entry[0]][2]       
            
            





            ##print "DAFUQ " + effallele + '\t' + dict_trait[entry[0]][1] + dict_trait[entry[0]][2] 
            #pass
   # pdb.set_trace()

    #risk allele alignment check
    



#pdb.set_trace()

#Go through dict and calculate wBs and w2s-2 for each snp, add up 
xtot = 0
ytot = 0
path = "./analysis"
if not os.path.exists(path):
    os.makedirs(path)
handle_indiv_result = file(trait_file + "_MR_result_indiv","w")
handle_indiv_result.write("SNP" + '\t' + "a_hat" + '\t' +"se_a" + '\t' + "chisq" +  '\t' + "p_value" + '\n')
for key in dict_trait.keys():
    if not '[NR]' in dict_trait[key]:
        if len(dict_trait[key]) == 5:
            templist = dict_trait[key]
            beta_trait = templist[0] 
            beta_disease = templist[3]
            se = templist[4]
                    ##print "seorig" + str(se)
            x = float(beta_trait) * float(beta_disease) * math.pow(float(se),-2)
               # #print "beta_trait" + str(float(beta_trait))
                ##print "beta_disease" +str(float(beta_disease))
               # #print "se" + str(math.pow(float(se),-2))
                #print x 
            y = math.pow(float(beta_trait),2) * math.pow(float(se),-2)
             #INDIV PART
            #pdb.set_trace()
            try:
                a_hat = float(x)/float(y)
                sea = float(math.sqrt(1/y))
                chisq = float(math.pow(a_hat/sea,2))
                p_value = 1 - stats.chi2.cdf(chisq,1)
                handle_indiv_result.write(key + '\t' + str(a_hat) + '\t' +str(sea) + '\t' + str(chisq) +  '\t' + str(p_value) + '\n') 
            except:
                print key + " is problematic"
                #ADD FIX FOR IF EFFECT IS 0, toss out earlier
                sys.exit()
            xtot += x
    	    ytot += y
        else:
            print "This SNP is not in the disease data : " + key
            
    
       
##print "xtot is " + str(xtot)
##print "ytot is " + str(ytot)

a_hat = float(xtot)/float(ytot)
sea = float(math.sqrt(1/ytot))
chisq = float(math.pow((a_hat/sea),2))
p_value = 1 - stats.chi2.cdf(chisq,1) #degrees of freedom always one?
#p_value figure out how to do scipy/numpy?

handle_result = file(trait_file + "_MR_result","w")
handle_result.write("a_hat = " + str(a_hat) + '\n' + "se(a) = " +str(sea) + '\n' + "chi_sq = " + str(chisq) + 
    '\n' + "p_value = " + str(p_value)) 

disease_handle.close()
trait_handle.close()
handle_result.close()


