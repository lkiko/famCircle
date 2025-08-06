# -*- encoding: utf-8 -*-
'''
@File        :line.py
@Time        :2021/07/11 12:44:00
@Author        :charles kiko
@Version        :1.0
@Contact        :charles_kiko@163.com
@Desc        :共线性局部
'''

import csv
import sys
import re
from math import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import *
from matplotlib.patches import Circle, Ellipse, Arc
from pylab import *
from famCircle.bez import *
plt.rc('font',family='Arial')

class line():
    def __init__(self, options):
        self.makers = 200
        self.class1 = False
        self.dpi = 600
        self.block = 6
        self.block_gep = 200
        self.position = 'order'
        self.GAP = 4
        self.h1,self.h2,self.h = 0.3,.6,.02
        self.blast_reverse = 'False'
        self.chr_name_size = '12'
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)
        self.block = int(self.block)
        self.block_gep = int(self.block_gep)
        self.GAP = float(self.GAP)
        self.h1 = float(self.h1)
        self.h2 = float(self.h2)
        self.h = float(self.h)


    def drow_coli(self,s1hpop,s2hpop,s1dpop,s2dpop):
        x,y = [],[]
        x.append(s1hpop)
        y.append(self.h1+(self.h)*1.1)
        x.append(s1dpop)
        y.append(self.h1+(self.h)*1.1)
        xt,yt = plot_bezier_curve(s1dpop,self.h1+(self.h)*1.1,s2dpop,self.h2-(self.h)*0.1)
        x = x + [i for i in xt]
        y = y + [i for i in yt]
        x.append(s2dpop)
        y.append(self.h2-(self.h)*0.1)
        x.append(s2hpop)
        y.append(self.h2-(self.h)*0.1)
        xt,yt = plot_bezier_curve(s2hpop,self.h2-(self.h)*0.1,s1hpop,self.h1+(self.h)*1.1)
        x = x + [i for i in xt]
        y = y + [i for i in yt]
        return x,y

    def run(self):
        lens1,chrlist1 = read_lens(self.lens1)
        lens2,chrlist2 = read_lens(self.lens2)
        # chrlist = chrlist1+chrlist2
        # lens = dict( lens1, **lens2)
        chrlist1,chrlist2 = self.chr1_name.split(','),self.chr2_name.split(',')
        gff1 = read_gff(self.gff1,chrlist1)
        gff2 = read_gff(self.gff2,chrlist2)
        gff = dict( gff1, **gff2)
        # 只保留需要的染色体共线性
        colineartly = readcolineartly(self.genepairs,self.block,self.block_gep,self.genepairsfile_type,gff,chrlist1+chrlist2)

        sumchr1,sumchr2 = sum([lens1[i][self.position] for i in chrlist1]),sum([lens2[i][self.position] for i in chrlist2])
        total_size = ((self.GAP+1)*max(sumchr1,sumchr2)/self.GAP)
        step = 0.99/total_size
        n1,n2 = len(chrlist1),len(chrlist2)
        se_dic1,se_dic2 = {},{}
        GAP1 = (total_size-sumchr1)/(n1+1)
        start = GAP1
        for j in range(n1):
            i = chrlist1[j]
            dic = {}
            dic['start'] = start
            dic['end'] = start + lens1[i][self.position]
            start = start + lens1[i][self.position] + GAP1
            se_dic1[i] = dic
            
        GAP2 = (total_size-sumchr2)/(n2+1)
        start = GAP2
        for j in range(n2):
            i = chrlist2[j]
            dic = {}
            dic['start'] = start
            dic['end'] = start + lens2[i][self.position]
            start = start + lens2[i][self.position] + GAP2
            se_dic2[i] = dic

        fig, ax = plt.subplots(figsize=(10, 10))
        ax.set_aspect('equal')
        
        for chro in chrlist1:
            xmin,xmax = se_dic1[chro]['start']*step,se_dic1[chro]['end']*step
            x,y = line_chr(self.h1,xmin,xmax,self.h)
            plt.text((xmin+xmax)/2, self.h1-0.02,chro,ha='center',va='center',fontsize=int(self.chr_name_size))
            plt.fill(x,y, facecolor=self.col1,alpha=.7)

        for chro in chrlist2:
            xmin,xmax = se_dic2[chro]['start']*step,se_dic2[chro]['end']*step
            x,y = line_chr(self.h2,xmin,xmax,self.h)
            plt.fill(x, y, facecolor=self.col2,alpha=.7)
            plt.text((xmin+xmax)/2, self.h2+self.h+0.02,chro,ha='center',va='center',fontsize=int(self.chr_name_size))
        for s1h,s2h,s1d,s2d in colineartly:
            if self.blast_reverse == "True":
                s1h,s2h = s2h,s1h
                s1d,s2d = s2d,s1d
            s1hpop,s2hpop,s1dpop,s2dpop = gff[s1h][self.position],gff[s2h][self.position],gff[s1d][self.position],gff[s2d][self.position]
            s1hchr,s2hchr,s1dchr,s2dchr = gff[s1h]['chr'],gff[s2h]['chr'],gff[s1d]['chr'],gff[s2d]['chr']
            if s1hchr not in (chrlist1+chrlist2) or s2hchr not in (chrlist1+chrlist2):# 排除不在绘制范围内的共线性
                continue
            if s1hchr != s2hchr:# 如果当前block在不同染色体上
                if s1hchr in chrlist1 and s2dchr in chrlist2:# 共线性第一列就是第一排
                    s1hpo,s2hpo,s1dpo,s2dpo = se_dic1[s1hchr]['start']+s1hpop,se_dic2[s2hchr]['start']+s2hpop,se_dic1[s1dchr]['start']+s1dpop,se_dic2[s2dchr]['start']+s2dpop
                    s1hpop,s2hpop,s1dpop,s2dpop = s1hpo*step,s2hpo*step,s1dpo*step,s2dpo*step
                    x,y = self.drow_coli(s1hpop,s2hpop,s1dpop,s2dpop)
                    plt.fill(x, y, facecolor=self.linec,alpha=.7)
                if self.genepairsfile_type == 'MCScanX':
                    if s1hchr in chrlist2 and s2dchr in chrlist1:# 共线性第一列就是第二排
                        s1h,s2h = s2h,s1h
                        s1d,s2d = s2d,s1d
                        s1hpop,s2hpop,s1dpop,s2dpop = gff[s1h][self.position],gff[s2h][self.position],gff[s1d][self.position],gff[s2d][self.position]
                        s1hchr,s2hchr,s1dchr,s2dchr = gff[s1h]['chr'],gff[s2h]['chr'],gff[s1d]['chr'],gff[s2d]['chr']         
                        s1hpo,s2hpo,s1dpo,s2dpo = se_dic1[s1hchr]['start']+s1hpop,se_dic2[s2hchr]['start']+s2hpop,se_dic1[s1dchr]['start']+s1dpop,se_dic2[s2dchr]['start']+s2dpop
                        s1hpop,s2hpop,s1dpop,s2dpop = s1hpo*step,s2hpo*step,s1dpo*step,s2dpo*step
                        x,y = self.drow_coli(s1hpop,s2hpop,s1dpop,s2dpop)
                        plt.fill(x, y, facecolor=self.linec,alpha=.7)

            elif s1hchr == s2hchr and s1hchr in chrlist1 and s2dchr in chrlist2:
                s1hpo,s2hpo,s1dpo,s2dpo = se_dic1[s1hchr]['start']+s1hpop,se_dic2[s2hchr]['start']+s2hpop,se_dic1[s1dchr]['start']+s1dpop,se_dic2[s2dchr]['start']+s2dpop
                s1hpop,s2hpop,s1dpop,s2dpop = s1hpo*step,s2hpo*step,s1dpo*step,s2dpo*step
                x,y = self.drow_coli(s1hpop,s2hpop,s1dpop,s2dpop)
                plt.fill(x, y, facecolor=self.linec,alpha=.7)
        # 关闭坐标轴
        ax.axis('off')
        plt.savefig(self.savefile,dpi=1000, bbox_inches = 'tight')