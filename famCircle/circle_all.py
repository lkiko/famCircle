# -*- encoding: utf-8 -*-
'''
@File        :circle_all.py
@Time        :2021/09/28 11:25:59
@Author        :charles kiko
@Version        :1.0
@Contact        :charles_kiko@163.com
@Desc        :基因关系圈图
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
from itertools import product
plt.rc('font',family='Arial')


class circle_all():
    def __init__(self, options):
    ##### circle parameters
        self.gap_ratio = 4 #gaps between chromosome circle, chr:gap = 4: 1
        self.radius = 0.3
        self.chrscize = 0.01
        self.dpi = 1000
        self.block_gep = 200
        self.genepairsfile_type = 'WGDI'
        self.chr_reverse = 'False'
        self.block = '0'
        self.class1 = True
        self.start_list = {}
        self.charc = 'None'
        self.color1 = 'None'
        self.seabornc = "rainbow"
        # self.seabornc = "cool"
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)
        self.gap_ratio = float(self.gap_ratio)
        self.chrscize = float(self.chrscize)

    def ksrun(self):
        lens1,chrlist1 = read_lens(self.lens1)
        lens2,chrlist2 = read_lens(self.lens2)
        cmap = make_cmap(self.seabornc)
        chro = chrlist1
        # if self.blast_reverse == 'True':
        #     chro = chrlist2
        colineartlyC = {}
        if self.color1 == 'None':
            for i in range(len(chro)):
                dic = {}
                dic['line'] = value_to_color(i, cmap, 0, len(chro))
                if self.charc == 'None':
                    dic['chr'] = value_to_color(i, cmap, 0, len(chro))
                else:
                    dic['chr'] = '#b3ada0'
                colineartlyC[chro[i]] = dic
        else:
            lt = self.color1.split(',')
            for i in lt:
                dic = {}
                dic['line'] = i.split(':')[2]
                if self.charc == 'None':
                    dic['chr'] = i.split(':')[1]
                else:
                    dic['chr'] = '#b3ada0'
                colineartlyC[i.split(':')[0]] = dic

        if chrlist1 == chrlist2:
            chrlist = chrlist1
        else:
            chrlist = chrlist1 + chrlist2


        lens = dict( lens1, **lens2)
        gff1 = read_gff(self.gff1,chrlist)
        gff2 = read_gff(self.gff2,chrlist)
        gff = dict( gff1, **gff2)
        colineartly = readcolineartly(self.genepairs,int(self.block),int(self.block_gep),self.genepairsfile_type,gff,chrlist)
        return gff,chrlist,lens,colineartly,colineartlyC

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


    def zj(self,lens,chr_list):
        fullchrolen = sum([lens[i]['end'] for i in lens.keys()])
        fullgene = sum([lens[i]['order'] for i in lens.keys()])
        gene_average = int(fullchrolen/fullgene)
        chr_number = len(lens.keys()) # total number of chromosomes
        GAP = fullchrolen/float(self.gap_ratio)/chr_number # gap size in base pair
        total_size = fullchrolen + chr_number * GAP # base pairs 
        for i in chr_list:
            # print(i)
            if i == chr_list[0]:
                self.start_list[i] = 0
            else:
                self.start_list[i] = self.start_list[chr_list[chr_list.index(i)-1]] + lens[chr_list[chr_list.index(i)-1]]['end'] + GAP
        return total_size,gene_average

    def run(self):
        # print("2023.8.8修订！")
        self.radius_a =  float(self.radius)
        self.radius_b = self.radius_a + float(self.chrscize) #.33, .335   # 半径控制参数 
        self.sm_radius=(self.radius_b-self.radius_a)/2 #telomere capping
        gff,chr_list,lens,colineartly,colineartlyC = self.ksrun()
        # print(chr_list)
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
            color = '#b3ada0'
            if i in colineartlyC.keys():
                color = colineartlyC[i]['chr']
            plt.fill(x, y, facecolor=color,alpha=.7)


            angle_deg = np.degrees((start+stop)/2)
            text_angle = angle_deg

            if 90 <= angle_deg <= 270:
                text_angle += 180
                alignment = 'right'
            else:
                alignment = 'left'


            label_x, label_y = self.rad_to_coord((start+stop)/2, self.radius_b*1.07)#1.2
            text(label_x, label_y, i, rotation=text_angle, horizontalalignment="center", verticalalignment="center", color = 'black', fontproperties=font_prop, fontsize = 7)

        up = []
        for local in colineartly:
            id1,id2,id3,id4 = local[0],local[1],local[2],local[3]

            if self.blast_reverse == 'True':
                id1,id2 = id2,id1
                id3,id4 = id4,id3

            chroa,chrob = gff[id1]['chr'],gff[id2]['chr']
            if chroa == chrob:
                # alp = 1
                up.append(local)
                continue
            else:
                alp = 0.5
            order1,order2,order3,order4 = gff[id1]['order'],gff[id2]['order'],gff[id3]['order'],gff[id4]['order']
            pos1,pos2,pos3,pos4 = gff[id1]['end'],gff[id2]['end'],gff[id3]['end'],gff[id4]['end']
            # col = self.color[chr_list.index(chroa)]

            color = '#b3ada0'
            if chroa in colineartlyC.keys():
                color = colineartlyC[chroa]['line']

            x,y = plot_colineartly_block(self.start_list,[chroa,chrob],[pos1,pos2,pos3,pos4], self.radius_a*0.99,total_size,self.chr_reverse,lens)
            plt.fill(x, y, facecolor=color,alpha=alp)
            # plot(x, y, '-', color='red', lw=2, alpha = .7)

        for local in up:
            # print(local)
            id1,id2,id3,id4 = local[0],local[1],local[2],local[3]
            chroa,chrob = gff[id1]['chr'],gff[id2]['chr']
            alp = 1
            order1,order2,order3,order4 = gff[id1]['order'],gff[id2]['order'],gff[id3]['order'],gff[id4]['order']
            pos1,pos2,pos3,pos4 = gff[id1]['end'],gff[id2]['end'],gff[id3]['end'],gff[id4]['end']
            color = '#b3ada0'
            if chroa in colineartlyC.keys():
                color = colineartlyC[chroa]['line']
            x,y = plot_colineartly_block(self.start_list,[chroa,chrob],[pos1,pos2,pos3,pos4], self.radius_a*0.99,total_size,self.chr_reverse,lens)
            # print(x,y,col)
            plt.fill(x, y, facecolor=color,alpha=alp)

            plot(x, y, '-.', color='black', lw=0.1, alpha = .3)

        plt.axis('off')
        savefig(self.savefile, dpi = int(self.dpi), bbox_inches = 'tight')# 
        sys.exit(0)
