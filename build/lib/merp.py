#!/bin/python

import requests
import pdb
import re
import os 
import sys
import math
from scipy import stats
#For plot visualization
import matplotlib.transforms as mtransforms
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
from pylab import *

class Merp():
    def pull(self,keyword=""):
        gwaswrite_handle = file('./data/finalgwas.txt', 'w')
        r = requests.get('http://www.genome.gov/admin/gwascatalog.txt')
        text = r.text
        unix_result = text.replace('\r', '')
        encoded_result = unix_result.encode('utf-8')
        gwaswrite_handle.write(encoded_result)
        gwaswrite_handle.close()
        gwasread_handle = file('./data/finalgwas.txt',"r")
        resultread_handle = file("result.txt","w")
        line = gwasread_handle.readlines()
        for entry in line[1:-1]: 
            entry_split = entry.split('\t')
            p = entry_split[27]
            '''Key Word Check'''
            if keyword in entry_split[7]:
                try:
                    p = float(p)
                    if p<=5.00E-8:
                        m = re.search('(?<=-)\w',entry_split[20])
                        allele = m.group()
                        if allele == "?":
                            allele = [NR]
                        resultread_handle.write(entry_split[13]+'\t'+entry_split[7]+'\t'+entry_split[27]+'\t'+entry_split[30]+ '\t'+ entry_split[31] +'\t'+entry_split[21]+'\t'+allele+'\t'+entry_split[26]+'\t'+entry_split[1]+'\n')
                except:
                    pass  



        resultread_handle.close()
        gwasread_handle.close()
        resultread_handle = file("result.txt","r")
        dict_trait = {}
        content = resultread_handle.readlines()
        _digits = re.compile('\d')
        _brackets = re.compile('\[.+?\]')
        var_dict = {}
        for line in content:
            entry = line.rstrip('\n').split('\t')
            if not(entry[1] in dict_trait):
                dict_trait[entry[1]] = []
            ci = entry[4]
            var_dict["ci_old"] = entry[4].strip()
            var_dict["ci_new"] = entry[4].strip()
            digits = bool(_digits.search(var_dict['ci_new']))
            brackets = bool(_brackets.search(var_dict['ci_new']))       
            #normalize cs so they all have brackets
            unit = '[NR]'
            '''Units Handling'''
            if digits:
                if 'NR' in var_dict['ci_new']:
                    var_dict['ci_new'] = '[NR]'
                if brackets:         
                    m = _brackets.match(var_dict['ci_new'])
                    brac = m.group(0)
                    var_dict['ci_new'] = brac.strip()
                    unit = var_dict["ci_old"].replace(var_dict['ci_new'],"")
                    unit = unit.strip()
                    #keep the ci column the same as before including units bc update will need it
                    ci = var_dict['ci_old']
            else:
                if ci == "NR" or ci == "NS":
                    ci = '[NR]'
                    unit = '[NR]'
                if '[NR]' not in ci:
                    if not ci == "":
                        ci = '[NR]'
                        unit = ci
                    else:
                        ci = '[NR]'
                        unit = '[NR]'
                if '[NR]' in ci:
                    unit = var_dict["ci_old"].replace('[NR]',"")
                    unit = unit.strip()
            #if no units given, fill with NR    
            if unit == '' or unit == "NR" or unit == "NS" :
                unit = '[NR]'
            if ci == '' or ci == "NR" or ci == "NS":
                ci = '[NR]'
            #Don't remove spaces in rs bc in filter will have to match up with nhgri otherwise key error
            dict_trait[entry[1]].append((entry[5],entry[3].replace(" ",""),unit,entry[6].replace(" ",""),entry[2].replace(" ",""),ci,entry[1],entry[0].replace(" ",""),entry[7].replace(" ",""),entry[8].replace(" ",""))) 
        handle_index = file("index","w")
        
        counter = 0
        path ='./traitFiles/'
        if not os.path.exists(path):
            os.makedirs(path)
        header = 'SNPrsID' +'\t' + 'Beta/OR' +'\t' +'Units' + '\t'+ 'Risk_Allele' +'\t' + 'p_val' +'\t' + '95%_CI' +'\t' + 'Trait'+ '\t' +'Gene'+'\t' +  'Risk_Allele_Freq' +'\t' +'PubMedID' +'\n'
        for key in dict_trait.keys():
            handle = file(path+str(counter),"w")
            templist = dict_trait[key]
            handle.write(header)
            snp_count = 0
            for tuple in templist:
                for string in tuple:
                	#Replace spaces and strip edges
                    comp_string = string.strip()
                    new_string = comp_string.replace (" ","_")
                    handle.write(str(new_string)+'\t')
                handle.write('\n')
                snp_count+=1
            handle.close()
            handle_index.write(str(counter)+'\t'+key+ '\t'+ str(snp_count) + " SNPs" + '\n')
            counter += 1
        handle_index.close()

    '''Helper for update'''
    def get_complement(self,nuc):
        if nuc == "A":
            return "T"
        if nuc == "T":
            return "A"
        if nuc == "G":
            return "C"
        if nuc == "C":
            return "G"

    def update_help(self,trait,allele_lines):
        handle = file(trait,"r")
        lines = handle.readlines()
        snp_dict = {}
        negative_list = ["decrease","lower","shorter","less"]
        write_handle = file(trait+"_update", "w")
        header = 'SNPrsID' +'\t' + 'Beta/OR' +'\t' +'Units' + '\t'+ 'Non-Risk_Allele' + '\t' 'Risk_Allele' +'\t' + 'p_val' +'\t' + '95%_CI' +'\t' + 'Trait'+ '\t' +'Gene'+'\t' +  'Risk_Allele_Freq' +'\t' +'PubMedID' +'\n'
        write_handle.write(header)
        for line in lines[1:]:
            entry = line.rstrip('\n').split('\t')
            risk_allele = entry[3]
            snp = entry[0]
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
                if risk == allele_one:
                    non_risk_dict[snp] = allele_two
                elif risk == allele_two:
                    non_risk_dict[snp] = allele_one
                else:
                    #if neither allele is the risk, take complements and work from there
                    complement = Merp.get_complement(self,risk)
                    if complement == allele_one:
                        non_risk_dict[snp] = Merp.get_complement(self,allele_two)
                    if complement == allele_two:
                        non_risk_dict[snp] = Merp.get_complement(self,allele_one)

        handle = file(trait,"r")
        lines = handle.readlines()
        not_in_genomes = ""
        for line in lines[1:]:
            entry = line.rstrip('\n').split('\t')
            risk_allele = entry[3]
            snp = entry[0]
            #Could also use units column, but CI contains both
            unit = entry[2]
            snp_dict[snp] = risk_allele
            try:
                allele_entry = non_risk_dict[snp]
            except:
                allele_entry = "[NR]" 
                not_in_genomes = not_in_genomes + snp + '\n'
            written = False
            for word in negative_list:
                if word in unit:
                    entry_unit = -(float(entry[1]))
                    beta_new = str(entry_unit)
                    write_handle.write(entry[0] + '\t' + beta_new + '\t' + entry[2] + '\t'+ allele_entry + '\t' + entry[3] + '\t' + entry[4]
                        + '\t' + entry[5]+ '\t'  + entry[6]+ '\t' + entry[7]+ '\t' + entry[8]+ '\t' + entry[9] + '\n' )
                    written = True
            if written == False:
                write_handle.write(entry[0] + '\t' + entry[1] + '\t' + entry[2] + '\t'+ allele_entry+ '\t' + entry[3] + '\t' + entry[4] + '\t' + entry[5] + '\t' + entry[6]+ '\t' + entry[7]+ '\t' + entry[8]+ '\t' + entry[9] + '\n' )
        if not_in_genomes != "":
            print not_in_genomes + "NON-RISK ALLELE NOT FOUND IN 1000GENOMES DATA FOR ABOVE SNPs FOR " + trait +  '.\n' +" PLEASE MANUALLY ENTER APPROPRIATE NUCLEOTIDE IN TEXT EDITOR USING OTHER RESOURCES OR DELETE SNP LINE " 


    def update(self,*traits):
        trait_list = []
        if len(traits) == 0:
            print "Usage: python update.py [trait]"
            sys.exit()
        for trait in traits:
            trait_list.append(trait)
        # allele_file = "data/1000_genomes"
        r = requests.get('http://coruscant.itmat.upenn.edu/merp/1000_genomes', stream=True)
        allele_lines = r.iter_lines()
        negative_list = ["decrease","lower","shorter"]
        nucleotides = ['A','T','G','C']
        # allele_handle = file(allele_file, "r")
        # allele_lines = allele_handle.readlines()
        allele_dict = {}
        for trait in trait_list:
            Merp.update_help(self,trait,allele_lines)
        ##loop back through trait file and add beta sign changing along with rewriting all lines 
   

    '''Helper functions for filter'''

    def trait_included(self,trait,include_list):
        for e in include_list:
            if e.lower() in trait.lower():
                return True
        return False

    def nhgri_test(self,snp,include_list,list_of_traits):
        mark = True
        for trait in list_of_traits:
            if not Merp.trait_included(self,trait,include_list):
                print snp + " has NHGRI catalog association with " + trait + " and is not exempt through nhgri_similar.txt. Throwing out."
                mark = False

        return mark

    def most_sig(self,snp_list,pval_dict):
        sig_snp = snp_list[0]
        for e in snp_list:
            if float(pval_dict[e]) <= float(pval_dict[sig_snp]):
                sig_snp = e
        return sig_snp

    #######pval violation >0.05 by chance algorithm#####
    # def dict_count(self,d):
    #     count = 0
    #     for key in d:
    #         #If passes conditions count it
    #         if d[key][0] == True:
    #             count += 1
    #     print count
    #     return count
    #Help function takes dictionary and cutoff and turns all keys with num_sig greater/= to cutoff to false
    #Then returns number of sig
    def dict_purge_count(self,d,cluster_d,cut,old_num_sig):
        #pdb.set_trace()
        return_sig = 0
        for key in d:
            if key in cluster_d:
                if d[key][1] >= cut:
                    cluster_d[key] = False 
            #get new number of sig
        ##currently ignore if not in pval file  
        for snp in d.keys():
            
            if snp in d: #unsure if dict_snp or dict
                if snp in cluster_d:
                    if cluster_d[snp] == True: 
                        return_sig = return_sig + d[snp][1]

        # for snp in snp_list:
        #   if snp == "rs4826508":
        #       pdb.set_trace()
        #   if key in d: #unsure if dict_snp or dict
        #       if key in cluster_d:
        #           if cluster_d[key] == True: 
        #               return_sig = return_sig + d[snp][1]

        print str(return_sig) + "is new sig after round changing cut"
        return return_sig

    def unit_checker(self,trait_file):
        #check this out, says unit in unit dict but not...
        unit_tags = ["decrease","lower","shorter","less","more","increase","higher","taller","greater"]
        unit_dict = {}
        handle = file(trait_file,"r")
        lines = handle.readlines()
        for entry in lines: 
            entry_split = entry.split('\t')
            unit = entry_split[2]
            if unit != "[NR]":
                unit = unit.lower()
            for tag in unit_tags:
                if tag in unit:
                    unit = unit.replace(tag,"*")
            if unit not in unit_dict.keys():
                unit_dict[unit] = []
        if len(unit_dict.keys()) > 1:
            print "File successsfuly filtered. Final filtered file can be found as traitFiles/" + trait_file +"filtered"
            print "WARNING: Different units detected in final filtered file. Please make desired unit conversions in a text editor:"
            for key in unit_dict.keys():
                print key
        else:
            print "File successfully filtered. Final filtered file can be found at" + trait_file +"filtered"
        return 0

    def filter(self,trait_file,include_file,exclude_file,out=False):
        if "update" not in trait_file:
            print "Warning: check to make sure you are using the updated file tagged with '_update' if using the pipeline"
            return
        if not Merp.file_checker(self,trait_file):
            return
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



        include_list = []
        exclude_list = []
        '''How to add usage error? '''
        #num_args = len(sys.argv) - 1
        # try:
        #   trait_file = sys.argv[1]
        #   include_file = sys.argv[2]
        #   exclude_file = sys.argv[3]
            
        # except:
        #   print "Usage: python filter.py [trait_file] [include_traits_file] [exclude_file]"
        #   sys.exit()
        #Read include list file and exclude list file(pval headers) 

        ###GETTING LD FILE FROM LD SNAP BROAD PROXY USING REQUESTS######
        rs_list = []
        with open(trait_file,"r") as trait:
            lines = trait.readlines()
            for line in lines[1:]:
                entry = line.rstrip('\n').split('\t')
                rs_list.append(entry[0])
        rs = ('\n\t').join(rs_list)
        rs =  '\n\t' +rs 

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
            snp = line_list[0]
            if snp not in snp_list:
                snp_list.append(snp)
        trait_handle.close()
        print str(len(snp_list)) + " is number of unique SNPs in the IVF entering the filter . . ."
        nhgri_handle = file(nhgri_file, "r")
        nhgri_lines = nhgri_handle.readlines()
        header = trait_lines[0]
        nhgri_dict = {}
        for line in nhgri_lines[1:]:
            #handles extra new line at end
            if line == '\n':
                break
            mod_line = line.rstrip('\n').split('\t')
            rs = mod_line[21]
            trait = mod_line[7]
            if rs in snp_list:
                if rs in nhgri_dict.keys() and trait not in nhgri_dict[rs]:
                    nhgri_dict[rs].append(trait)
                else:
                    nhgri_dict[rs] = [trait]
        nhgri_handle.close()


        snp_nhgri_dict = {}
        for snp in snp_list:
            if snp in nhgri_dict.keys():
                list_of_traits = nhgri_dict[snp]
                snp_nhgri_dict[snp] = Merp.nhgri_test(self,snp,include_list,list_of_traits)  
            elif out==False:
                print snp + " is not found in NHGRI catalog. If you want you are bringing outside SNPs into filter, please add out=True as argument"
                snp_nhgri_dict[snp] = False
                pass
            elif out == True:
                print snp + " is not in NHGRI catalog but kept in because out=True. Please be cautious for confounding."



        print '\n'
        ####NHGRI PORTION END#######
                
        ######PVAL PORTION START########
        # pval_handle = file(pval_file, "r")
        p = requests.get('http://coruscant.itmat.upenn.edu/merp/v3abr_allmetabolic_pvals_v2.txt', stream=True)
        pval_lines = p.iter_lines()
        dict_snp = {}
        # pval_lines = pval_handle.readlines()
        #header = pval_lines[0]
        # header =  "SNP CHR POS P_SBP P_DBP P_HDL P_LDL P_TG P_BMI P_HEIGHT P_WHRadjBMI P_2hrGLUadjBMI P_FastGlu P_HbA1C P_FastIns P_HOMA-B P_HOMA-IR P_FastProIns"
        header = pval_lines.next()
        header_split = header.rstrip('\n').split(' ')
        columns = len(header_split) - 3
       # correction_num = 0
        excluded_index = []
        for p in header_split:
            for excluded_trait in exclude_list:
                if excluded_trait.lower() in p.lower():
                    print p + " associations from our metabolic confounding filter will be ignored because found in pval_similar.txt"
                    #will add the index of the column we want to exclude to a list 
                    excluded_index.append(header_split.index(p))
                    columns = columns - 1
                    #correction_num += 1

        '''Review why removed replacement pval thingy'''            
        # replace_handle = file(replacement_pval_snps, 'r')
        # replace_lines = replace_handle.readlines()
        # replace_dict = {}
        # for line in replace_lines[1:]:
        #   entry = line.rstrip('\n').split('\t')
        #   orig = entry[0]
        #   replacement = entry[1]
        #   if orig in no_pval_snps.keys():
        #       if orig not in replace_dict.keys():
        #           replace_dict[orig] = replacement
        #       else:
        #           print orig + " is repeated in replacement text file"
        #           pass     
        for line in pval_lines:
            mod_line = line.rstrip('\n').split(' ')
            rs = mod_line[0]
            if rs in snp_list:
                count_p1 = 0
                count_p2 = 0
                for e in mod_line[3:]:
                    if mod_line.index(e) in excluded_index:
                        continue
                    if float(e) <= pmax1:
                        count_p1 = count_p1 + 1
                    if float(e) <= pmax2:
                        count_p2 = count_p2 + 1
            # count_p1 = count_p1 - correction_num
            # count_p2 = count_p2 - correction_num
            #only add to dict if rs in trait file
                if rs not in dict_snp.keys(): 
                    #first compare true count to modified threshold then change count to modified count if trait_related is true
                    if count_p1 <= threshold1:
                        if count_p2 <= threshold2:
                            #if trait in pval file, subtract 1 from count of number of sig associations
                            dict_snp[rs] = [True, count_p1] # Record number of associations wiht p<0.05 afte rfiltering for less than 4
                        else:
                            dict_snp[rs] = [False, count_p1]
                    else:
                        dict_snp[rs] = [False, count_p1]
            
            '''Replacement pval code'''
            # elif rs in replace_dict.keys():
            #   if count_p1 <= threshold1:
            #       if count_p2 <= threshold2:
            #           #if trait in pval file, subtract 1 from count of number of sig associations
            #           dict_snp[rs] = [True, count_p1] # Record number of associations wiht p<0.05 afte rfiltering for less than 4
            #       else:
            #           dict_snp[rs] = [False, count_p1]
            #   else:
            #       dict_snp[rs] = [False, count_p1]
        # pval_handle.close()

        '''IF NOT IN PVAL### keep in snp_list for now, later change so replaces?'''

        no_pval_snps = {}
        for snp in snp_list:
            if snp not in dict_snp.keys():
                no_pval_snps[snp] = []
        #pdb.set_trace()
        #Calculate total number of significant assoc (<0.05)
        num_sig = 0
        trait_handle = file(trait_file,"r")
        trait_lines = trait_handle.readlines()
        for line in trait_lines[1:]:
            line_list = line.rstrip('\n').split('\t')
            if len(line_list) > 5:
                snp = line_list[0]
            else:
                print "file format error."
            if snp in dict_snp:
                #Create new dict with only trait files
                if dict_snp[snp][0] == True: 
                    num_sig = num_sig + dict_snp[snp][1]
        trait_handle.close()

        ''' PVAL PORTION END '''

        '''LD PORTION BEGIN'''

        trait_handle = file(trait_file,"r")
        pval_dict = {}
        repeated_snps = []
        trait_lines = trait_handle.readlines()
        for line in trait_lines[1:]:
            line_list = line.rstrip('\n').split('\t')
            snp = line_list[0]
            p = line_list[5]
            if snp not in pval_dict:
                pval_dict[snp] = p
            else:
                '''Repeating SNPs handling'''

                p_new = p
                #replace p only if more sig?
                if float(p_new) < float(pval_dict[snp]):
                    pval_dict[snp] = p_new
                if snp not in repeated_snps:
                    repeated_snps.append(snp)
        trait_handle.close()

        #Go through and rewrite traitfile without repeating snps

        '''Repeating SNPs handling'''
        trait_handle = file(trait_file,"r")
        trait_lines = trait_handle.readlines()
        header = trait_lines[0]
        abridged_trait_handle = file("./"+trait_file + "_abr_temp","w")
        abridged_trait_handle.write(header)
        for line in trait_lines[1:]:
            line_list = line.split('\t')
            snp = line_list[0]
            p = line_list[5]
            pid = line_list[10]
            if snp not in repeated_snps:
                abridged_trait_handle.write(line)
            else:
                if float(p) == float(pval_dict[snp]):
                    abridged_trait_handle.write(line)

        abridged_trait_handle.close()
        trait_handle.close()

        cluster_dict = {}
        index_dict = {}

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

        '''Review necessity of this'''   
        if len(not_in_ld)>=1:
            for line in ld_result_hm22[1:]:
                line_list = line.split('\t')
                if len(line_list) != 1:
                    snp = line_list[0]
                    proxy = line_list[1]
                    if snp in not_in_ld:
                        if "WARNING" in proxy: #and "query" in proxy #and "not" in proxy and "in" in proxy:
                            if "Query snp not in" in line_list[2]: 
                                if snp not in not_in_ld:
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

        for snp in not_in_ld:
            for key in index_dict.keys():
                if snp==key: #test
                    not_in_ld.remove(snp)

        abr_trait_handle = file(trait_file+"_abr_temp","r")
        trait_lines = abr_trait_handle.readlines()
        header = trait_lines[0]
        cluster_status = {}
        most_sig_cluster = {}

        for line in trait_lines[1:]:
            line_list = line.split('\t')
            snp = line_list[0]
            p = line_list[5]
            if snp in not_in_ld:
                pass
            elif snp in index_dict.keys():
                index = index_dict[snp]
                assoc_list = cluster_dict.get(index)
                sig_rs = Merp.most_sig(self,assoc_list,pval_dict)

                most_sig_cluster[index] = sig_rs
                '''Check this'''
                #Case where most sig snp not found in pval file
                # case = False
                # pdb.set_trace()
                # while case == False:
                #   if sig_rs not in dict_snp.keys():
                #       sig_rs = most_sig(assoc_list.remove(sig_rs))
                #   if sig_rs in dict_snp.keys():
                #       case == True

                # if case == False:
                #   pass





                '''eventually delete'''
                print sig_rs + " is the most sig out of cluster " + str(index)
        abr_trait_handle.close()
        for key in most_sig_cluster.keys():
            snp = most_sig_cluster[key]
            index = index_dict[snp]
            assoc_list = cluster_dict.get(index)
            assoc_list.remove(snp)
            #nhgri test
            if snp_nhgri_dict[snp] == True:
                ##pval test
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
                        #   rsq = ld_assoc_dict[(snp,rs)]

                        #   if rsq >= rsq_nhgri_threshold:
                        #       if snp_nhgri_dict[rs] == False:
                        #           cluster_status[snp] = False
                        #           break
                        #       else: 
                        #           pass
                        #           #elif rsq >= rsq_nhgri_lower_threshold:
                        #           #   if snp_nhgri_dict[rs] == False:
                        #           #       cluster_status[snp] = "IDK"
                        #           #       print "ambigious for" + snp
                        #           #       pass
                        #           #   else: 
                        #           cluster_status[snp] = True
                        #   else:
                        #       cluster_status[snp] = True
                    else:
                        cluster_status[snp] = False
                else:
                    cluster_status[snp] = False
            else:
                cluster_status[snp] = False
            ''' If snp is not in pval file, keep in and put warning'''
            if snp in no_pval_snps.keys():
                print snp + ': WARNING: SNP not found in in pval_assoc.txt.' + '\n' + 'Keeping in for now- please double-check associations with this SNP before proceeding with analysis' 

                cluster_status[snp] = True

        threshold_met = False
        new_num_sig = num_sig
        cut = threshold1
        while threshold_met == False:
            ####JUST LOOK AT ELIGIBLE SNPS
            num_snps = 0
            for key in cluster_status.keys():
                if cluster_status[key] == True:
                    num_snps +=1
            num_tests = num_snps * columns
            if float(new_num_sig) <= max_fraction * num_tests:
                threshold_met = True
            if float(new_num_sig) > max_fraction * num_tests:
                old_num_sig = new_num_sig
                new_num_sig = Merp.dict_purge_count(self,dict_snp,cluster_status,cut,old_num_sig)
                cut = cut - 1
                if cut < 0:
                    print "Fatal error: cut below 0, consider having a higher cutoff threshold of violations or higher max_fraction "
                    sys.exit()
        snps_to_write = []
        non_assoc_write_temp = []
        for key in cluster_status.keys():
            if cluster_status[key] == True:
                snps_to_write.append(key)
        '''Try to get path to work, keep trying'''
        # path="./traitFiles/final_trait_files"
        # if not os.path.exists(path):
        #     os.makedirs(path)
        abr_trait_handle = file(trait_file+"_abr_temp","r")
        trait_lines = abr_trait_handle.readlines()
        header = trait_lines[0]
        updated_handle = file(trait_file +"filtered", "w")
        updated_handle.write(header)
        for line in trait_lines[1:]:
            line_list = line.split('\t')
            snp = line_list[0]
            if snp in snps_to_write:
                updated_handle.write(line)
                '''Do we really not care about SNPs not in LD data?'''
            if snp in not_in_ld:
                print snp + ' is not in LD data search for 1000genomes or Hm22. We reccommend tossing this SNP. Alternatively, you may manually check for LD'
                updated_handle.write('---Not in LD Data---' + line)



        os.remove(trait_file+"_abr_temp")
        abr_trait_handle.close()

        Merp.unit_checker(self,trait_file)

        #If A and B are both in dict 1 and significant assoc, change b index to A's in dict 2, add B and all clusters with B
        # to index of A's cluster list and remove index of B's entry of cluster lists.
        #If A or B is in dict 1 but the other is not, simply add the new snp to the cluster index in dict 1 and index in dict 2
        #If neither A or B are in dict 1, create a new cluster and index in dict 1 and ref in dict 2

    def hasNumbers(self,instring):
        return any(char.isdigit() for char in instring)

    def file_checker(self,trait_file):
        nuc = ['A','T','G','C']
        with open(trait_file,"r") as infile:
            lines = infile.readlines()
            header = lines[0]
            header = header.split('\t')
            if "rs" not in header[0].lower() or "beta" not in header[1].lower() or "unit" not in header[2].lower() or "non-risk_allele" not in header[3].lower() or "risk_allele" not in header[4].lower():
                print "Header incorrectly formatted. Please make sure your file is formatted in the following manner:" + '\n' + "SNPrsID   Beta/OR Units   Non-Risk_Allele    Risk_Allele . . . ."
                return
            for line in lines[1:]:
                entry = line.split('\t')
                if len(entry) <=1:
                    print "Extra new lines detected. Please remove excess new lines."
                    return False
                rs = entry[0]
                beta = entry[1]
                unit = entry[2]
                non_risk = entry[3]
                risk = entry[4]
                if "rs" not in rs:
                    print rs + " is not a valid SNPrsID. Please edit"
                    return False
                if not Merp.hasNumbers(self,beta):
                    print beta + " is not a valid beta value. Please enter a numerical value for SNP " + rs 
                    return False
                if non_risk not in nuc:
                    print non_risk + " is not a valid Non-Risk Allele entry. Please use A, T, G or C"
                    return False
                if risk not in nuc:
                    print risk  +" is not a valid Non-Risk Allele entry. Please use A, T, G or C"
                    return False
            print "File format check PASSED!"
            return True






    def calc(self,trait_file,disease_file):
        
        

        # try:
        #     trait_file= sys.argv[1]
        #     disease_file = sys.argv[2]
        # except:
        #     print "Usage: python v3mr_calc_fin.py [path/to/traitfile] [path/to/diseasefile]"
        #     sys.exit()

        ######TRAIT FILE PARAMETERS######
        rs_index = 0
        beta_index = 1
        unit_index = 2
        allele1_index = 3
        riskallele_index = 4


        ####DISEASE FILE PARAMTERS#######
        dis_rs_index = 0
        dis_allele1_index = 1
        dis_riskallele_index = 2
        dis_lnOR_index = 3
        dis_lnse_index = 4

        if not Merp.file_checker(self,trait_file):
            return
        print 'Calculating effect now . . .'


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
        # path = "./analysis"
        # if not os.path.exists(path):
        #     os.makedirs(path)
        handle_indiv_result = file(trait_file + "_MR_result_indiv","w")
        handle_indiv_result.write("SNP" + '\t' + "a_hat" + '\t' +"se_a" + '\t' + "chisq" +  '\t' + "p_value" + '\n')
        for key in dict_trait.keys():
            if not '[NR]' in dict_trait[key]:
                if len(dict_trait[key]) == 5:
                    templist = dict_trait[key]
                    beta_trait = templist[0] 
                    beta_disease = templist[3]
                    se = templist[4]
                    #For T2D format
                    # OR_disease = templist[3]
                    # beta_disease = math.log(float(OR_disease))
                    # LCI = templist[4]
                    # se_temp = float(math.log(float(OR_disease))) - float(math.log(float(LCI)))
                    # se = se_temp/1.96
                    x = float(beta_trait) * float(beta_disease) * math.pow(float(se),-2)

                    y = math.pow(float(beta_trait),2) * math.pow(float(se),-2)
                     #INDIV PART
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
        print "ADD SOME LOG SUMMARY + WHERE FILES ARE"

    def plot(self,indiv_file,summ_file,Trait,Disease):
        

        # try:
        #   indiv_file = sys.argv[1]
        #   summ_file = sys.argv[2]
        #   Trait = sys.argv[3]
        #   Disease = sys.argv[4]
            
        # except:
        #   print "Usage: python calc_to_plot.py [indiv_file] [summ_file] [Trait] [Disease]"
        #   sys.exit()
        indiv_handle = file(indiv_file, "r")
        summ_handle = file(summ_file, "r")

        indiv_lines = indiv_handle.readlines()
        indiv_dict = {}
        for line in indiv_lines[1:]:
            entry = line.rstrip('\n').split('\t')
            snp = entry[0]
            a_hat = entry[1]
            se_a = entry[2]
            indiv_dict[snp] = [a_hat,se_a]

        plot_list = []
        y = 2
        for key in indiv_dict:
            #pdb.set_trace()

            a_hat = indiv_dict[key][0]
            se = indiv_dict[key][1]
            l95 = float(a_hat) - 1.96*float(se)
            new_se = float(a_hat) - l95

            temp_list = [key,a_hat, y, new_se]
            y += 1
            plot_list.append(temp_list)

        figsize = (6,6)
        fig = figure(figsize=figsize, dpi=80)
        ax_kms = SubplotHost(fig, 1,1,1, aspect=1.)

        fig.add_subplot(ax_kms)

        ## Add overall OR first 
        summ_lines = summ_handle.readlines()
        summ_list =[]
        for line in summ_lines:
            entry = line.rstrip('\n').split(' ')
            summ_list.append(entry[2])
        #pdb.set_trace()
        ov_a_hat = summ_list[0]
        ov_se = summ_list[1]
        # ov_odds   = math.exp(float(ov_a_hat))

        ##
        ov_l95 = float(ov_a_hat) - 1.96*float(ov_se)
        ov_new_se = float(ov_a_hat) - ov_l95

        ax_kms.errorbar(float(ov_a_hat), 1, xerr=ov_new_se, color = "r", fmt = 'o')

        ##
        for key, a, y, new_se in plot_list:
            ax_kms.errorbar([float(a)], y, xerr=[new_se], color="k", fmt='o')


        ax_kms.axis["bottom"].set_label("Log of Odds Ratio")
        ax_kms.axis["left"].set_label("SNPs")
        title('Estimated Causal Effect of ' + Trait + ' on ' + Disease)
        rank_list = [1]
        rs_list = ['Overall']
        count = 2  
        for snp in indiv_dict.keys():
            rs_list.append(snp)
            rank_list.append(count)
            count +=1
        #if Trait == "Height":
        #   rs_list = ['Overall']

        plt.yticks(rank_list,rs_list)
        plt.xticks([-2,-1,0,1,2],[-2,-1,0,1,2] )
        ax_kms.set_ylim(0, len(rs_list) + 1)
        axvline(x=0)


        show()

# if __name__ == '__main__':
#     #invoke the freeze_support funtion for windows based systems
#     try:
#         sys.getwindowsversion()
#         multiprocessing.freeze_support()
#     except:
#         pass

#     x = Cluster()
#     exit()



