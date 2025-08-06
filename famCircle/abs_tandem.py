#absolute（绝对）、relative（相对）
import numpy as np
import pandas as pd
from tqdm import trange
from famCircle.bez import *
import matplotlib.pyplot as plt

class abs_tandem():
	def __init__(self,options):
		self.gap = 30
		# self.Orthogroups = 'Orthogroups.tsv'
		# self.Orthogroups_UnassignedGenes = 'Orthogroups_UnassignedGenes.tsv'
		# self.gff = 'sir_gene.gff'
		# self.tandem_file = 'sir.tandem'
		for k, v in options:
			setattr(self, str(k), v)
			print(k, ' = ', v)

	def hist_(self,ax,list0,name,color,binx):
	    ax.grid(c='grey',ls='--',linewidth=0.3)
	    # nums = np.array(list0)
	    # list0 = list(nums[nums<=0.4])
	    # nums = np.array(list0)
	    # list0 = list(nums[nums>0])
	    ax.hist(list0,binx,alpha=0.4,color=color,label=name,edgecolor=color)
	    ax.tick_params(labelsize=15)
	    ax.legend(fontsize=15)

	def read_gff(self,gff):
		gene = {}
		gff = pd.read_csv(gff,header = None, sep='\t', comment='#')
		for index,row in gff.iterrows():
			dic0 = {}
			dic0['chr'] = row[0]
			dic0['end'] = row[2]
			dic0['order'] = row[5]
			gene[row[1]] = dic0
		return gene

	def read_Orthogroups(self):
		Orthogroups = pd.read_csv(self.orthogroups,header = 0, sep='\t', comment='#')
		print(Orthogroups)
		tandem = []
		tandem_gap = []
		tandem_gapend = []
		gff = self.read_gff(self.gff)
		for index,row in Orthogroups.iterrows():
			if row[self.genome] is not np.nan:
				print(row['Orthogroup'])
				lt = row[self.genome].split(', ')
				lt_chr = [gff[i]['chr'] for i in lt]
				dic = {}
				for i in trange(len(lt_chr)):
					if lt_chr[i] not in dic.keys():
						dic[lt_chr[i]] = [lt[i]]
					else:
						dic[lt_chr[i]].append(lt[i])
				print(dic)
				for i in dic.keys():
					lt = dic[i]
					lt_index = [gff[x]['order'] for x in lt]
					np_matx = np.zeros((len(lt),len(lt)))
					print(np_matx)
					for x in range(len(lt)):
						for y in range(x+1,len(lt)):
							np_matx[x][y] = abs(lt_index[x] - lt_index[y])
					print(np_matx)
					
					for x in range(len(lt)):
						for y in range(x+1,len(lt)):
							if np_matx[x][y] <= int(self.gap):
								print('tandem ',lt[x],lt[y])
								tandem.append(lt[x])
								tandem.append(lt[y])
								tandem_gap.append(abs(gff[lt[x]]['order']-gff[lt[y]]['order']))
								tandem_gapend.append(abs(gff[lt[x]]['end']-gff[lt[y]]['end']))

		print(len(list(set(tandem))))
		tandem = list(set(tandem))
		f=open(self.tandem_file,'w')
		f.write('\n'.join(tandem))
		f.close()

		fig,ax=plt.subplots(2,1,figsize=(20,15),dpi=800)
		self.hist_(ax[0],tandem_gap,'tandem_order','#38b48b',int(self.gap))
		self.hist_(ax[1],tandem_gapend,'tandem_end','#83ccd2',100)

		# plt.legend(fontsize=15)
		plt.xlabel('gap',fontsize=20)# 设置 x 轴标签
		plt.ylabel('density',fontsize=20)# 设置 y 轴标签

		plt.savefig(self.savefile, dpi = 1000, bbox_inches = 'tight')# 存储图片

	def read_Orthogroups_UnassignedGenes(self):
		Orthogroups_UnassignedGenes = pd.read_csv(self.orthogroups_UnassignedGenes,header = 0, sep='\t', comment='#')
		print(Orthogroups_UnassignedGenes)

	def run(self):
		self.read_Orthogroups()

# a = abs_tandem()
# a.read_Orthogroups()
# # a.read_Orthogroups_UnassignedGenes()