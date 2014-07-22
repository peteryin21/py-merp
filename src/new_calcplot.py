#!/bin/python
import pdb
import matplotlib.transforms as mtransforms
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
from pylab import *
import sys
try:
	indiv_file = sys.argv[1]
	summ_file = sys.argv[2]
	Trait = sys.argv[3]
	Disease = sys.argv[4]
	
except:
	print "Usage: python calc_to_plot.py [indiv_file] [summ_file] [Trait] [Disease]"
	sys.exit()
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
# ov_odds	= math.exp(float(ov_a_hat))

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
#	rs_list = ['Overall']

plt.yticks(rank_list,rs_list)
plt.xticks([-2,-1,0,1,2],[-2,-1,0,1,2] )
ax_kms.set_ylim(0, len(rs_list) + 1)
axvline(x=0)


show()
