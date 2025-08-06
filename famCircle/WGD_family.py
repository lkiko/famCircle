# -*- encoding: utf-8 -*-
'''
@File        :tamdem.py
@Time        :2021/09/28 11:17:18
@Author        :charles kiko
@Version        :1.0
@Contact        :charles_kiko@163.com
@Desc        :放射型圈图
'''

import re
import os
import sys
import gc
from math import *
import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import *
from pylab import *
from famCircle.bez import *
from matplotlib import gridspec
import itertools
plt.rc('font',family='Arial')


class WGD_family():
    def __init__(self, options):
    ##### circle parameters
        self.gap_ratio = 5 #gaps between chromosome circle, chr:gap = 4: 1
        self.radius = 0.3
        self.dpi = 600
        self.block= 0.008 #block scize
        self.blockthick = 0.007 #0.006
        self.start_list = {}
        self.class1 = True
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)

    def ksrun(self):
        lens,chr_list = read_lens(self.lens)# 染色体字典
        gff = read_gff(self.gff)# 基因字典
        blast = self.readblast(gff)
        ks = read_ks0(self.ks)# ks字典
        family = read_family(self.genefamily)# 家族列表

        # for i in gff.keys():
        #     if gff[i]['chr'] not in lens.keys():
        #         gff = gff.pop(i)
        for i in family:# 去除不在染色体上的家族成员
            if gff[i]['chr'] not in lens.keys():
                family.remove(i)
        return lens,chr_list,gff,blast,ks,family


    def readblast(self,gff):
        one_gene = []
        alphagenepairs = open(self.genepairs, 'r', encoding='utf-8')
        if self.file_type == 'famCircle':
            for row in alphagenepairs:
                if (row[0] == '\n'):
                    continue
                elif (row[0] == '#'):
                    continue
                else:
                    lt = row.strip('\n').split()
                    id1, id2 = str(lt[0]),str(lt[1])
                    if (id1 not in gff.keys() or id2 not in gff.keys()):
                        continue
                    one_gene.append(id1)
                    one_gene.append(id2)
        elif self.file_type == 'WGDI':
            for row in alphagenepairs:
                if (row[0] == '\n'):
                    continue
                elif (row[0] == '#'):
                    continue
                else:
                    lt = row.strip('\n').split()
                    id1, id2 = str(lt[0]),str(lt[2])
                    if (id1 not in gff.keys() or id2 not in gff.keys()):
                        continue
                    one_gene.append(id1)
                    one_gene.append(id2)
        elif self.file_type == 'ColinearScan':
            for row in alphagenepairs:
                if (row[0] == '\n' and row[0] == '+'):
                    continue
                elif (row[:3] == 'the'):
                    continue
                else:
                    lt = row.strip('\n').split()
                    id1, id2 = str(lt[0]),str(lt[2])
                    if (id1 not in gff.keys() or id2 not in gff.keys()):
                        continue
                    one_gene.append(id1)
                    one_gene.append(id2)
        elif self.file_type == 'MCScanX':
            for row in alphagenepairs:
                if (row[0] == '\n'):
                    continue
                elif(row[:12] == '## Alignment'):
                    continue
                elif ('#' not in row):
                    lt = row.strip('\n').split()
                    # print(lt)
                    id1, id2 = str(lt[2]),str(lt[3])
                    if (id1 not in gff.keys() or id2 not in gff.keys()):
                        continue
                    one_gene.append(id1)
                    one_gene.append(id2)
        else:
            print('genepairsfile_type error: File Format not recognized!')
            exit()
        # self.colornum = len(chrnum)
        alphagenepairs.close()
        # blastlist = sorted(one_gene,key=(lambda x:x[0]))
        # # print(one_gene)
        # blast = {}
        # for i in one_gene:
        #     if i[0] not in blast.keys():
        #         lt = []
        #         lt.append(i[1])
        #         blast[i[0]] = lt
        #     else:
        #         blast[i[0]].append(i[1])

        return one_gene

    def rad_to_coord(self, angle, radius):
        return radius*cos(angle), radius*sin(angle)

    def to_deg(self, bp, total):
        # from basepair return as degree
        return bp*360./total

    def to_radian(self, bp, total):
        # from basepair return as radian
        # print ("to_radian", bp, total)
        return radians(bp*360./total)

    def plot_arc(self, start, stop, radius):
        # start, stop measured in radian 染色体
        #print start, stop
        t = arange(start, stop, pi/720.)
        x, y = radius*cos(t), radius*sin(t)
        plot(x, y, "k-",alpha=.5)

    def plot_cap(self, angle, clockwise):
        radius=self.sm_radius
        # angle measured in radian, clockwise is boolean 鸭舌
        if clockwise: 
            t = arange(angle, angle+pi, pi/30.)
        else: 
            t = arange(angle, angle-pi, -pi/30.)
        x, y = radius*cos(t), radius*sin(t)
        middle_r = (self.radius_a+self.radius_b)/2
        x, y = x + middle_r*cos(angle), y + middle_r*sin(angle)
        plot(x, y, "k-", alpha=.5)

    def plot_arc_block(self, start, radius):# block
        t = arange(start, start+self.block, pi/720.)
        # print(start, start+self.block)
        x,y = radius * cos(t), radius*sin(t)
        x1, y1 = (radius-self.blockthick) * cos(t), (radius-self.blockthick) * sin(t)
        plot(x, y, "b-", linewidth=3, alpha=0.9)

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

    def transform_deg(self, ch, pos, total_size):
        return self.to_deg(pos + self.start_list[ch], total_size)

    def transform_pt(self, ch, pos, r, total_size):
        # convert chromosome position to axis coords
        rad = self.to_radian(pos + self.start_list[ch], total_size)
        return r*cos(rad), r*sin(rad)

    def transform_pt1(self,ch, pos, r, total_size):
        # convert chromosome position to axis coords
        # print("transform", ch, pos, r)
    #    print "startlist", self.start_list[ch]
        rad = self.to_radian(pos + self.start_list[ch], total_size)
        return r*cos(rad), r*sin(rad),rad

    def transform_pt2(self, rad, r):
        return r*cos(rad), r*sin(rad)

    def plot_bez_inner(self, p1, p2, cl, total_size,lw0):
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
        xt = Bezier(x, t)# 贝塞尔曲线
        yt = Bezier(y, t)
        plot(xt, yt, '-', color=cl, lw=lw0, alpha=0.3)#alpha 透明度

    def plot_bez_Ks2(self, rad1, r1, rad2, r2, col, ratio):
    #    print "bez Ks 2"
        # print('rad1, r1',rad1, r1)
        ex1x, ex1y = self.transform_pt2(rad1, r1)
        ex2x, ex2y = self.transform_pt2(rad2, r2)
        # ratio = -0.7#0.5
        sita = pi / 2
        if ex1x != ex2x:
            sita = atan((ex2y-ex1y)/(ex2x-ex1x))
        d = sqrt((ex2x-ex1x)**2+(ex2y-ex1y)**2)
        L = d * ratio
        P1x = ex1x + L*sin(sita)
        P1y = ex1y - L*cos(sita)
        P2x = ex2x + L*sin(sita)
        P2y = ex2y - L*cos(sita)
        step = .01
        t = arange(0, 1+step, step)
        x=[ex1x, P1x, P2x, ex2x]
        y=[ex1y, P1y, P2y, ex2y]
        # print('x,y,t',x,y,t)
        xt = Bezier(x,t)
        yt = Bezier(y,t)
        plot(xt, yt, '-', color = col, lw = 0.7)#0.1

    def run(self):
        self.radius_a =  float(self.radius)
        self.radius_b = self.radius_a + 0.005#.33, .335   # 半径控制参数 
        self.sm_radius=(self.radius_b-self.radius_a)/2 #telomere capping
        if (os.path.exists(self.outfile)):
            os.remove(self.outfile)
        lens,chr_list,gff,blast,ks,family = self.ksrun()
        total_size,gene_average = self.zj(lens,chr_list)
        stop_list = {}

        fig = plt.figure(figsize=(10, 11))
        gs = gridspec.GridSpec(30, 30)#, width_ratios=[10, 1]
        # gs0 = gridspec.GridSpec(1, 2)
        ax1 = fig.add_subplot(gs[1:25,1:25], aspect=1)

        for i in self.start_list.keys():
            stop_list[i] = self.start_list[i] + lens[i]['end']
        for i in chr_list:
            start,stop = self.start_list[i],stop_list[i]
            start, stop = self.to_radian(start, total_size), self.to_radian(stop, total_size)
            # shaft
            self.plot_arc(start, stop, self.radius_a)
            self.plot_arc(start, stop, self.radius_b)
            # telemere capping
            clockwise=False
            self.plot_cap(start, clockwise)
            clockwise=True
            # print(start)
            self.plot_cap(stop, clockwise)
            label_x, label_y = self.rad_to_coord((start+stop)/2, self.radius_b*1.07)#1.2
            #print label_x, label_y
            text(label_x, label_y, i, horizontalalignment="center", verticalalignment="center", fontsize = 7, color = 'black')

        # 绘制关系，并查询WGD_family
        pair = list(itertools.permutations(family,2))
        tp_list = []
        for i in pair:
            id1,id2 = i[0],i[1]
            if id1 == id2:
                continue
            chro1,chro2 = gff[id1]['chr'],gff[id2]['chr']
            if chro1 == chro2 and abs(int(gff[id1]['order'])-int(gff[id2]['order'])) == 1:
                if id1 not in tp_list:
                    tp_list.append(id1)
                if id2 not in tp_list:
                    tp_list.append(id2)
        tandem_list = tp_list
        prox_list = []
        for i in pair:
            id1,id2 = i[0],i[1]
            if id1 == id2:
                continue
            pos1,pos2 = gff[id1]['end'],gff[id2]['end']
            chro1,chro2 = gff[id1]['chr'],gff[id2]['chr']
            if chro1 == chro2 and id1 not in tandem_list and id2 in tandem_list and abs(int(gff[id1]['order'])-int(gff[id2]['order'])) < 15:
                if id1 not in prox_list:
                    prox_list.append(id1)
            elif chro1 == chro2 and id1 in tandem_list and id2 not in tandem_list and abs(int(gff[id1]['order'])-int(gff[id2]['order'])) < 15:
                if id2 not in prox_list:
                    prox_list.append(id2)
        tp_list = tandem_list + prox_list
        outfile = open(self.outfile, 'w', encoding='utf-8')
        out_list = []
        for i in family:
            if i in blast and i not in tp_list:
                if i not in out_list:
                    outfile.write(str(i) + '\n')
                    out_list.append(i)
                    chro1 = gff[i]['chr']
                    if chro1 not in chr_list:
                        continue
                    posi = gff[i]['end'] + self.start_list[gff[i]['chr']]
                    start = self.to_radian(posi, total_size)
                    self.plot_arc_block(start, self.radius_b)# 基因
        outfile.close()
        del pair
        gc.collect()

        pair = list(itertools.permutations(out_list,2))

        pu_ = []
        for i in pair:
            id1,id2 = i[0],i[1]
            if id1 == id2:
                continue
            pos1,pos2 = gff[id1]['end'],gff[id2]['end']
            chro1,chro2 = gff[id1]['chr'],gff[id2]['chr']
            if str(id1)+'_'+str(id2) in ks.keys():
                ks_v = ks[str(id1)+'_'+str(id2)]
            elif str(id2)+'_'+str(id1) in ks.keys():
                ks_v = ks[str(id2)+'_'+str(id1)]
            else:
                continue
            if ks_v <= float(self.ks_concern.split(',')[0]) or ks_v >= float(self.ks_concern.split(',')[1]):
                self.plot_bez_inner((chro1, pos1, self.radius_a*0.95), (chro2, pos2, self.radius_a*0.95), "#e9dfe5", total_size,2)#绘制内部共线性
            else:
                pu_.append([id1,id2])
                continue
        for i in pu_:
            id1,id2 = i[0],i[1]
            pos1,pos2 = gff[id1]['end'],gff[id2]['end']
            chro1,chro2 = gff[id1]['chr'],gff[id2]['chr']
            if str(id1)+'_'+str(id2) in ks.keys():
                ks_v = ks[str(id1)+'_'+str(id2)]
                col = return_col(ks_v,float(self.ks_concern.split(',')[0]),float(self.ks_concern.split(',')[1]))
            elif str(id2)+'_'+str(id1) in ks.keys():
                ks_v = ks[str(id2)+'_'+str(id1)]
                col = return_col(ks_v,float(self.ks_concern.split(',')[0]),float(self.ks_concern.split(',')[1]))
            else:
                continue
            self.plot_bez_inner((chro1, pos1, self.radius_a*0.98), (chro2, pos2, self.radius_a*0.98), col, total_size,2)#绘制内部共线性
        plt.axis('off')
        drawSpec(float(self.ks_concern.split(',')[0]),float(self.ks_concern.split(',')[1]),gs[:,25:])
        savefig(self.savefile, dpi = int(self.dpi), bbox_inches = 'tight')
        sys.exit(0)
