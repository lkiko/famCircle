# -*- encoding: utf-8 -*-
'''
@File		:circle_all.py
@Time		:2021/09/28 11:25:59
@Author		:charles kiko
@Version		:1.0
@Contact		:charles_kiko@163.com
@Desc		:基因关系圈图
'''


import re
import sys
from math import *
import gc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import *
from matplotlib.patches import Circle, Ellipse
from pylab import *
from collections import Counter
from famCircle.bez import *
from famCircle.bez import drw_line, get_genedensity, get_genomeGC, plot_colineartly_block, rad_to_coord, read_fasta, readcolineartly,plot_arc_local

plt.rc('font',family='Arial')


class re_loops():
	def __init__(self, options):
	##### circle parameters
		self.gap_ratio = 2 #gaps between chromosome circle, chr:gap = 4: 1
		self.gap_r = 0.95 #右侧信息窗口占比
		self.radius = 0.5
		self.chrscize = 0.007
		self.dpi = 1000
		self.block_gep = 200
		self.genepairsfile_type = 'MCScanX'
		self.draw_type = "Coverage"
		self.block = '0'
		self.class1 = True
		self.chr_reverse = "False"
		self.loopnsize = 0.03
		self.start_list = {}
		for k, v in options:
			setattr(self, str(k), v)
			print(k, ' = ', v)

	def ksrun(self):
		lens1,chrlist = read_lens(self.lens)
		reflens,refchrlist = read_lens(self.reflens)
		gff1 = read_gff(self.gff,chrlist)
		refgff = read_gff(self.refgff,refchrlist)
		gff = dict( gff1, **refgff)
		lens = dict( lens1, **reflens)

		colineartly = readcolineartly(self.genepairs,int(self.block),int(self.block_gep),self.genepairsfile_type,gff,chrlist+refchrlist)
		repeat,repeatsof = [],[]
		gene_Coverage,gene_overlap = get_genedensity(self.gff,lens1,int(self.windows),int(self.step),chrlist,2,3)
		# print(gene_Coverage)
		GC = {}
		if self.genomefile != 'None':
			genome = read_fasta(self.genomefile)
			GC = get_genomeGC(genome,int(self.windows),int(self.step),chrlist)
		
		if self.repeatfile != "None":
			lt = self.repeatfile.split(',')
			repeatsof = self.repeatsof.split(',')
			for i in range(len(lt)):
				Coverage,overlap = get_genedensity(lt[i],lens1,int(self.windows),int(self.step),chrlist,3,4)
				repeat.append(Coverage)
				repeat.append(overlap)
		repeat.append(gene_Coverage)
		repeat.append(gene_overlap)
		repeatsof.append('GENE')
		return gff,chrlist,lens1,lens,refchrlist,reflens,colineartly,GC,repeatsof,repeat

	def rad_to_coord(self, angle, radius):
		return radius*cos(angle), radius*sin(angle)

	def to_radian(self, bp, total):
		# from basepair return as radian
		return radians(bp*360./total)

	def plot_arc(self, start, stop, radius):
		# start, stop measured in radian
		t = arange(start, stop, pi/720.)
		x, y = radius*cos(t), radius*sin(t)
		plot(x, y, "k-", alpha=.5)# 染色体圆弧

	def plot_cap(self, angle, clockwise):
		radius=self.sm_radius
		# angle measured in radian, clockwise is boolean
		if clockwise: 
			t = arange(angle, angle+pi, pi/30.)
		else: 
			t = arange(angle, angle-pi, -pi/30.)
		x, y = radius*cos(t), radius*sin(t)
		middle_r = (self.radius_a+self.radius_b)/2
		x, y = x + middle_r*cos(angle), y + middle_r*sin(angle)
		plot(x, y, "k-", alpha=.5)# 边缘

	def zj(self,lens,chr_list,reflens,refchr_list):
		fullchrolen = sum([lens[i]['end'] for i in lens.keys()])
		fullgene = sum([lens[i]['order'] for i in lens.keys()])
		gene_average = int(fullchrolen/fullgene)
		chr_number = len(lens.keys()) # total number of chromosomes
		GAP = fullchrolen/float(self.gap_ratio)/chr_number # gap size in base pair
		total_size = fullchrolen + chr_number * GAP # base pairs 
		GAP = GAP * float(self.gap_r)
		allchr = refchr_list + chr_list
		for i in allchr:
			# print(i)
			if i == allchr[0]:
				self.start_list[i] = 0
			else:
				self.start_list[i] = self.start_list[allchr[allchr.index(i)-1]] + lens[allchr[allchr.index(i)-1]]['end'] + GAP
		return total_size,gene_average,GAP

	def transform_pt(self, ch, pos, r, total_size):
		rad = self.to_radian(pos + self.start_list[ch], total_size)
		return r*cos(rad), r*sin(rad)

	def plot_bez_inner(self, p1, p2, cl, total_size,alp,lw):
	#	print "inner"
		a, b, c = p1
		# print(a,'1')
		ex1x, ex1y = self.transform_pt(a, b, c, total_size)
		a, b, c = p2
		# print(a,'2')
		ex2x, ex2y = self.transform_pt(a, b, c, total_size)
		# Bezier ratio, controls curve, lower ratio => closer to center
		ratio = .5
		x = [ex1x, ex1x*ratio, ex2x*ratio, ex2x]
		y = [ex1y, ex1y*ratio, ex2y*ratio, ex2y]
		step = .01
		t = arange(0, 1+step, step)
		xt = Bezier(x, t)
		yt = Bezier(y, t)
		plot(xt, yt, '-', color=cl, lw=lw, alpha = alp)#alpha 


	def draw_colineartly(self,coldic,gff,chr_list,lens,colineartly,total_size,gene_average,s_e_dic1):

		up = []
		for local in colineartly:

			id1,id2,id3,id4 = local[0],local[1],local[2],local[3]
			chroa,chrob = gff[id1]['chr'],gff[id2]['chr']
			if chroa == chrob:
				up.append(local)
				continue
			alp = 0.5
			pos1,pos2,pos3,pos4 = gff[id1]['end'],gff[id2]['end'],gff[id3]['end'],gff[id4]['end']
			# col = self.color[chr_list.index(chroa)]
			col = coldic[chrob]
			x,y = plot_colineartly_block(self.start_list,[chroa,chrob],[pos1,pos2,pos3,pos4], self.radius*0.99,total_size,self.chr_reverse,lens)
			plt.fill(x, y, facecolor=col,alpha=alp)
			# plot(x, y, '-', color='red', lw=2, alpha = .7)

		for local in up:
			# print(local)
			id1,id2,id3,id4 = local[0],local[1],local[2],local[3]
			chroa,chrob = gff[id1]['chr'],gff[id2]['chr']
			alp = 1
			# order1,order2,order3,order4 = gff[id1]['order'],gff[id2]['order'],gff[id3]['order'],gff[id4]['order']
			pos1,pos2,pos3,pos4 = gff[id1]['end'],gff[id2]['end'],gff[id3]['end'],gff[id4]['end']
			# col = self.color[chr_list.index(chroa)]
			col = coldic[chrob]
			x,y = plot_colineartly_block(self.start_list,[chroa,chrob],[pos1,pos2,pos3,pos4], self.radius*0.99,total_size,self.chr_reverse,lens)
			# print(x,y,col)
			plt.fill(x, y, facecolor=col,alpha=alp)
			plot(x, y, '-.', color='black', lw=0.5, alpha = .5)



	def read_data(self,soft,pvuv,x_index,y_index,chr_list,s_e_dic,gff,chr_dic,total_size,r_min,r_max,col,lens):
		size = float(self.loopnsize)/12
		drw_line(chr_list,s_e_dic,r_min,False,None,"b-",0.2,0.1)
		# drw_line(chr_list,s_e_dic,r_max,False,None,"b-",0.2,0.1)
		r_min,r_max = r_min + size,r_max - size
		data = pvuv.set_index([1])[2].to_dict()



		# avg = pvuv[y_index].mean()
		# print(pvuv)
		# max_d = pvuv[y_index].max()
		# min_d = pvuv[y_index].min()
		# avg_line = ((avg-min_d)/(max_d-min_d))*(r_max-r_min)+r_min
		# drw_line(chr_list,s_e_dic,avg_line,False,None,"g-",0.4,0.1)
		dic = pvuv.groupby(0).groups# 按照第几列分组
		for key in dic.keys():
			if key not in chr_list:
				continue
			print(" ******** ",key," ******** ")
			if key == chr_list[0]:
				label_x, label_y = rad_to_coord(-0.04, (r_min+r_max)/2)
				text(label_x, label_y, soft, ha='left',va='center',horizontalalignment="center", verticalalignment="center", fontsize = 3, color = 'black',rotation=-90)		 #ha='left',#x=2.2是文字的左端位置，可选'center', 'right', 'left' va='baseline',#y=8是文字的低端位置，可选'center', 'top', 'bottom', 'baseline', 'center_baseline'
			else:
				pass
			x1 = []# 低于均值
			y1 = []
			x2 = []# 高于均值
			y2 = []
			x3_1 = []# 均值线
			y3_1 = []
			x3_2 = []# 均值线
			y3_2 = []
			local = pvuv.loc[dic[key]].sort_values(by=[0,2],ascending= [True,True])

			avg = local[y_index].mean()
			max_d = local[y_index].max()
			min_d = local[y_index].min()
			avg_line = ((avg-min_d)/(max_d-min_d))*(r_max-r_min)+r_min

			starto = s_e_dic[key]['start']
			endo = s_e_dic[key]['end']
			plot_arc_local(starto, endo, avg_line,"g-",0.4,0.1)
			label_x, label_y = rad_to_coord(s_e_dic[key]['start']-0.01, avg_line)
			text(label_x, label_y, str(100*avg)[:5], horizontalalignment="center", verticalalignment="center", fontsize = 2, color = 'black',rotation=(s_e_dic[key]['start']/(2*pi))*360)



			for index,row in local.iterrows():
				if self.chr_reverse == 'False':
					start = data[row[x_index]]
				else:
					if lens[row[0]]['rev'] == '-':
						start = lens[row[0]]['end'] - data[row[x_index]]
					else:
						start = data[row[x_index]]


				x = s_e_dic[row[0]]['start'] + radians((start)*360.0/total_size)
				# x = s_e_dic[row[0]]['start'] + radians((data[row[x_index]])*360.0/total_size)
				# x = s_e_dic[row[0]]['start'] + radians((gff[row[x_index]])*360.0/total_size)
				y = ((row[y_index]-min_d)/(max_d-min_d))*(r_max-r_min)+r_min
				label_x, label_y = rad_to_coord(x, y)
				x1.append(label_x)
				y1.append(label_y)

			plt.plot(x1, y1, '-', color=col, alpha=1, linewidth=0.2, label='low')
			resx = x1
			for i in range(len(x3_1)-1,-1,-1):
				resx.append(x3_1[i])
			resy = y1
			for i in range(len(y3_1)-1,-1,-1):
				resy.append(y3_1[i])
			# plt.fill(resx ,resy,facecolor='green',alpha=0.5)






	def draw_loops(self,gff,chr_list,lens,GC,repeatsof,repeat,loopN,s_e_dic,total_size):
		cols = ['#2ca9e1','#eb6238','#68be8d','#867ba9','#b3ada0','#028760','#82ae46','#f39800','#8491c3','#bb5548']
		if self.genomefile != 'None':
			print('GC')
			self.read_data("GC",GC,1,7,chr_list,s_e_dic,gff,lens,total_size,self.radius-((len(repeatsof)+1)*float(self.loopnsize)),self.radius-((len(repeatsof))*float(self.loopnsize)),cols[0],lens)
		if self.draw_type == 'Coverage':
			for i in range(len(repeatsof)):
				print(repeatsof[i])
				self.read_data(repeatsof[i],repeat[i*2],1,7,chr_list,s_e_dic,gff,lens,total_size,self.radius-((i+1)*float(self.loopnsize)),self.radius-((i)*float(self.loopnsize)),cols[i+1],lens)
		elif self.draw_type == 'overlap':
			for i in range(len(repeatsof)):
				print(repeatsof[i])
				self.read_data(repeatsof[i],repeat[i*2+1],1,7,chr_list,s_e_dic,gff,lens,total_size,self.radius-((i+1)*float(self.loopnsize)),self.radius-((i)*float(self.loopnsize)),cols[i+1],lens)



	def run(self):
		# print("2023.8.8修订！")
		coldic = {}
		for i in self.refchr2col.split(','):
			coldic[i.split(':')[0]] = i.split(':')[1]

		gff,chr_list,lens1,lens,refchr_list,reflens,colineartly,GC,repeatsof,repeat = self.ksrun()
		loopN = 0
		if len(GC) != 0:
			loopN += 1
		if len(repeatsof) != 0:
			loopN = loopN + len(repeatsof)
		self.radius_a =  float(self.radius) - (loopN * float(self.loopnsize))
		# print(self.radius_a)
		self.radius_b = self.radius_a + float(self.chrscize) #.33, .335   # 半径控制参数 
		self.sm_radius=(self.radius_b-self.radius_a)/2 #telomere capping
		total_size,gene_average,gap_length = self.zj(lens,chr_list,reflens,refchr_list)

		s_e_dic = {}
		s_e_dic1 = {}
		start = 0.0
		allchr = refchr_list + chr_list
		for i in allchr:
			dic_0 = {}# 弧度制
			dic_0['start'] = radians(start*360.0/total_size)
			dic_0['end'] = radians((start + lens[i]['end'])*360.0/total_size)
			s_e_dic[i] = dic_0

			dic_01 = {}
			dic_01['start'] = int(start)
			dic_01['end'] = int(start + lens[i]['end'])
			s_e_dic1[i] = dic_01

			start = start + lens[i]['end'] + gap_length

		fig = plt.figure(figsize=(10, 10))
		for key in allchr:
			label_x, label_y = self.rad_to_coord((s_e_dic[key]['start']+s_e_dic[key]['end'])/2, self.radius*1.07)#1.2
			text(label_x, label_y, key, horizontalalignment="center", verticalalignment="center", fontsize = 7, color = 'black')
			start,stop = s_e_dic1[key]['start'],s_e_dic1[key]['end']
			drw_line(allchr,s_e_dic,0.51,False,None,"k-",0.2,1)
			# start, stop = self.to_radian(start, total_size), self.to_radian(stop, total_size)
			# # shaft
			start, stop = self.to_radian(start, total_size), self.to_radian(stop, total_size)
			x,y = plot_rec(start, stop,0.46,self.chrscize)
			if key in coldic.keys():
				col = coldic[key]
			else:
				col = 'black'
			plt.fill(x, y, facecolor=col,alpha=.7)

		self.draw_colineartly(coldic,gff,chr_list,lens,colineartly,total_size,gene_average,s_e_dic1)
		# self.draw_loops(gff,chr_list,lens,GC,repeatsof,repeat,loopN,s_e_dic,total_size)
		plt.axis('off')
		plt.savefig(self.savefile, dpi = int(self.dpi), bbox_inches = 'tight')
		# plt.show()
		sys.exit(0)