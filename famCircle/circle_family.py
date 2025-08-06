# -*- encoding: utf-8 -*-
'''
@File        :circle_family.py
@Time        :2021/09/28 11:25:41
@Author        :charles kiko
@Version        :1.0
@Contact        :charles_kiko@163.com
@Desc        :家族分布圈图
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
plt.rc('font',family='Arial')

class circle_family():
    def __init__(self, options):
    ##### circle parameters
        self.gap_ratio = 2 #gaps between chromosome circle, chr:gap = 4: 1
        self.radius = 0.3
        self.chrscize = 0.01
        self.dpi = 600
        self.genepairsfile_type = 'WGDI'
        self.block = '0'

        self.class1 = True
        self.start_list = {}
        self.color = ["#2ca9e1","#ffd900","#e2041b","#38b48b","#ea5506","#8f2e14"
                      ,"#b44c97","#f5b1aa","#bce2e8","#82ae46","#824880","#5a544b","#4d4c61"
                      ,"#6f4b3e","#192f60","#2b2b2b","#475950","#727171","#2c4f54","#455765","#432f2f"]
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)
        self.chrC,self.famC,self.bacC = self.chr_fam_bacc.split(',')[0],self.chr_fam_bacc.split(',')[1],self.chr_fam_bacc.split(',')[2]
        self.gap_ratio = float(self.gap_ratio)


    def ksrun(self):
        lens,chrlist = read_lens(self.lens)
        gff = read_gff(self.gff,chrlist)
        colineartly = self.readblast()
        family,family_dic,dic_motif = read_family1(self.genefamily)# 家族列表
        family = drop_none(family,gff)
        return gff,chrlist,lens,colineartly,family

    def rad_to_coord(self, angle, radius):
        return radius*cos(angle), radius*sin(angle)

    def to_radian(self, bp, total):
        # from basepair return as radian
        return radians(bp*360./total)

    def plot_arc(self, start, stop, radius):
        # start, stop measured in radian
        t = arange(start, stop, pi/720.)
        x, y = radius*cos(t), radius*sin(t)
        plot(x, y, '-', color=self.chrC, alpha=.5)# 染色体圆弧

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
        plot(x, y, '-', color=self.chrC, alpha=.5)# 边缘

    def zj(self,lens,chr_list):
        fullchrolen = sum([lens[i]['end'] for i in lens.keys()])
        fullgene = sum([lens[i]['order'] for i in lens.keys()])
        gene_average = int(fullchrolen/fullgene)
        chr_number = len(lens.keys()) # total number of chromosomes
        GAP = fullchrolen/float(self.gap_ratio)/chr_number # gap size in base pair
        total_size = fullchrolen + chr_number * GAP # base pairs 
        for i in chr_list:
            if i == chr_list[0]:
                self.start_list[i] = 0
            else:
                self.start_list[i] = self.start_list[chr_list[chr_list.index(i)-1]] + lens[chr_list[chr_list.index(i)-1]]['end'] + GAP
        return total_size,gene_average

    def transform_pt(self, ch, pos, r, total_size):
        rad = self.to_radian(pos + self.start_list[ch], total_size)
        return r*cos(rad), r*sin(rad)

    def transform_pt0(self, ch, pos, r, total_size):
        rad = self.to_radian(pos + self.start_list[ch], total_size)
        return r*cos(rad), r*sin(rad),rad

    def readblast(self):
        blast_dic = {}
        one_gene = []
        alphagenepairs = open(self.genepairs, 'r', encoding='utf-8')
        if self.genepairsfile_type == 'famCircle':
            for row in alphagenepairs:
                if (row[0] == '\n'):
                    continue
                elif (row[0] == '#'):
                    lt = row.strip('\n').split(' ')
                    block = {}
                    for i in lt:
                        if '=' in str(i):
                            lt0 = i.split('=')
                            block[str(lt0[0])] = str(lt0[1])
                    N = int(block['N'])
                    if N >= int(self.block):
                        self.class1 = True
                    else :
                        self.class1 = False
                else:
                    if self.class1:
                        lt = row.strip('\n').split()
                        id1, id2 = str(lt[0]),str(lt[1])
                        if id1 not in one_gene:
                            one_gene.append(id1)
                            blast_dic[id1] = [id2]
                        else:
                            blast_dic[id1].append(id2)
        elif self.genepairsfile_type == 'WGDI':
            for row in alphagenepairs:
                if (row[0] == '\n'):
                    continue
                elif (row[0] == '#'):
                    lt = row.strip('\n').split(' ')
                    block = {}
                    for i in lt:
                        if '=' in str(i):
                            lt0 = i.split('=')
                            block[str(lt0[0])] = str(lt0[1])
                    N = int(block['N'])
                    if N >= int(self.block):
                        self.class1 = True
                    else :
                        self.class1 = False
                else:
                    if self.class1:
                        lt = row.strip('\n').split()
                        id1, id2 = str(lt[0]),str(lt[2])
                        if id1 not in one_gene:
                            one_gene.append(id1)
                            blast_dic[id1] = [id2]
                        else:
                            blast_dic[id1].append(id2)
        elif self.genepairsfile_type == 'ColinearScan':
            for row in alphagenepairs:
                if (row[0] == '\n',row[0] == '+'):
                    continue
                elif (row[:3] == 'the'):
                    lt = row.strip('\n').split()
                    N = int(lt[-1])
                    if N >= int(self.block):
                        self.class1 = True
                    else :
                        self.class1 = False
                else:
                    if self.class1:
                        lt = row.strip('\n').split()
                        id1, id2 = str(lt[0]),str(lt[2])
                        if id1 not in one_gene:
                            one_gene.append(id1)
                            blast_dic[id1] = [id2]
                        else:
                            blast_dic[id1].append(id2)
        elif self.genepairsfile_type == 'MCScanX':
            for row in alphagenepairs:
                if (row[0] == '\n'):
                    continue
                elif(row[:12] == '## Alignment'):
                    lt = row.strip('\n').split(' ')
                    block = {}
                    for i in lt:
                        if '=' in str(i):
                            lt0 = i.split('=')
                            block[str(lt0[0])] = str(lt0[1])
                    N = int(block['N'])
                    if N >= int(self.block):
                        self.class1 = True
                    else :
                        self.class1 = False
                elif ('#' not in row):
                    if self.class1:
                        lt = row.strip('\n').split()
                        # print(lt)
                        if len(lt) == 5:
                            id1, id2 = str(lt[2]),str(lt[3])
                        elif len(lt) == 4:
                            id1, id2 = str(lt[1]),str(lt[2])
                        elif len(lt) == 6:
                            id1, id2 = str(lt[3]),str(lt[4])
                        else:
                            print(row)
                            print('Parse error!')
                            exit()
                        if id1 not in one_gene:
                            one_gene.append(id1)
                            blast_dic[id1] = [id2]
                        else:
                            blast_dic[id1].append(id2)
        elif self.genepairsfile_type == 'BLAST':
            num = 0
            name0_list = []
            for row in alphagenepairs:
                if (row[0] == '\n'):
                    continue
                elif ('#' not in row):
                    lt = row.strip('\n').split()
                    if lt[0] == lt[1]:
                        continue
                    if str(lt[0]) not in name0_list:
                        name0_list.append(str(lt[0]))
                        num = 1
                    else:
                        if num <= int(self.block):
                            pass
                        else:
                            continue
                    id1, id2 = str(lt[0]),str(lt[1])
                    if id1 not in one_gene:
                        one_gene.append(id1)
                        blast_dic[id1] = [id2]
                    else:
                        blast_dic[id1].append(id2)
        else:
            print('genepairsfile_type error: File Format not recognized!')
            exit()
        return blast_dic

    def plot_bez_inner(self, p1, p2, cl, total_size,alp,lw):
    #    print "inner"
        a, b, c = p1
        ex1x, ex1y = self.transform_pt(a, b, c, total_size)
        a, b, c = p2
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

    def run(self):
        self.radius_a =  float(self.radius)
        self.radius_b = self.radius_a + float(self.chrscize) #.33, .335   # 半径控制参数 
        self.sm_radius=(self.radius_b-self.radius_a)/2 #telomere capping
        gff,chr_list,lens,colineartly,family = self.ksrun()
        total_size,gene_average = self.zj(lens,chr_list)
        stop_list = {}
        for i in self.start_list.keys():
            stop_list[i] = self.start_list[i] + lens[i]['end']
        fig = plt.figure(figsize=(10, 10))
        path = get_path()
        font_path = os.path.join(path, 'example/arial.ttf')
        from matplotlib.font_manager import FontProperties
        font_prop = FontProperties(fname=font_path)
        for i in chr_list:
            start,stop = self.start_list[i],stop_list[i]
            start, stop = self.to_radian(start, total_size), self.to_radian(stop, total_size)
            # shaft
            x,y = plot_arc(start, stop,self.radius_a,self.chrscize)
            plt.fill(x, y, facecolor=self.chrC,alpha=.7)
            label_x, label_y = self.rad_to_coord((start+stop)/2, self.radius_b*1.07)#1.2
            text(label_x, label_y, i, horizontalalignment="center", verticalalignment="center", color = 'black', fontproperties=font_prop, fontsize = 7)
        famlist = []
        list0 = []
        dic0 = {}
        for id1 in colineartly.keys():
            if id1 not in gff.keys() or gff[id1]['chr'] not in lens.keys():
                continue
            for id2 in colineartly[id1]:
                if id2 not in gff.keys() or gff[id1]['chr'] not in lens.keys():
                    continue
                chro1,chro2 = gff[id1]['chr'],gff[id2]['chr']
                order1,order2 = gff[id1]['order'],gff[id2]['order']
                alp = 0.03
                if abs(order1-order2)<50:
                    pass
                else:
                    if id1 in family or id2 in family:
                        famlist.append([id1,id2])
                    else:
                        if chro1 not in list0:
                            list0.append(chro1)
                            dic0[chro1] = [[id1, id2]]
                        else:
                            dic0[chro1].append([id1, id2])
        for chro in dic0:
            lt = dic0[chro]
            for i in lt:
                id1, id2 = i[0],i[1]
                pos1,pos2 = gff[id1]['end'],gff[id2]['end']
                chro1,chro2 = gff[id1]['chr'],gff[id2]['chr']
                alp = 0.08
                self.plot_bez_inner((chro1, pos1, self.radius_a*0.8), (chro2, pos2, self.radius_a*0.8), self.bacC, total_size,alp,2)
        ltx = []
        for i in famlist:
            id1, id2 = i[0],i[1]
            pos1,pos2 = gff[id1]['end'],gff[id2]['end']
            chro1,chro2 = gff[id1]['chr'],gff[id2]['chr']
            alp = 1
            self.plot_bez_inner((chro1, pos1, self.radius_a*0.8), (chro2, pos2, self.radius_a*0.8), self.famC, total_size,alp,2.5)
            a, b, c = chro1, pos1, self.radius_a*0.91
            ex1x, ex1y,rad1 = self.transform_pt0(a, b, c, total_size)
            a, b, c = chro2, pos2, self.radius_a*0.91
            ex2x, ex2y,rad2 = self.transform_pt0(a, b, c, total_size)
            if id1 not in ltx:
                ltx.append(id1)
                # ex1y = ex1y + 0.004
                text(ex1x, ex1y, id1, horizontalalignment="center", verticalalignment="center", color = 'black',rotation=(rad1/(2*pi))*360, fontproperties=font_prop, fontsize = 6)
            if id2 not in ltx:
                ltx.append(id2)
                text(ex2x, ex2y, id2, horizontalalignment="center", verticalalignment="center", color = 'black',rotation=(rad2/(2*pi))*360, fontproperties=font_prop, fontsize = 6)
        plt.axis('off')
        savefig(self.savefile, dpi = int(self.dpi), bbox_inches = 'tight')
        sys.exit(0)
