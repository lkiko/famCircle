# -*- encoding: utf-8 -*-
'''
@File        :lineblock.py
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

class lineblock():
    def __init__(self, options):
        self.makers = 200
        self.class1 = False
        self.dpi = 600
        self.block = 5
        self.start = None
        self.end = None
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)

    def map1(self,axes1,maker,x,y,x1,color):
        # makers = 15 #标尺
        # x = 0.1 #起始坐标x
        # y = 0.1 #起始坐标y
        # x1 = 0.8#染色体长度
        # color='k'
        lw = 1  #比例线宽
        y1 = 0.001
        w = 0.01
        alpha = 1
        # print(y+w+y1, x, x+x1,y, x, x+x1, lw)
        plt.axhline(y=y, xmin=x, xmax=x+x1, lw=lw, c=color, alpha=alpha)
        plt.axhline(y=y+w+y1, xmin=x, xmax=x+x1, lw=lw, c=color, alpha=alpha)
        for i in range(maker+1):
            mx = x + ((x1/maker)*i)
            plt.axvline(x=mx, ymin=y, ymax=y+(w/2.5), lw=lw, c=color, alpha=alpha)
        base1 = Arc(xy=(x, y+((w+y1)/2)),    # 椭圆中心，（圆弧是椭圆的一部分而已）
                width=w+y1,    # 长半轴
                height=w+y1,    # 短半轴
                angle=90,    # 椭圆旋转角度（逆时针） 
                theta1=0,    # 圆弧的起点处角度
                theta2=180,    # 圆度的终点处角度
                color=color,
                alpha=alpha,
                linewidth=lw
                )
        base2 = Arc(xy=(x1+x, y+((w+y1)/2)),    # 椭圆中心，（圆弧是椭圆的一部分而已）
                width=w+y1,    # 长半轴
                height=w+y1,    # 短半轴
                angle=-90,    # 椭圆旋转角度（逆时针） 
                theta1=0,    # 圆弧的起点处角度
                theta2=180,    # 圆度的终点处角度
                color=color,
                alpha=alpha,
                linewidth=lw   #线宽像素
                )
        axes1.add_patch(base1)
        axes1.add_patch(base2)

    def map2(self,axes1,x,y,x1,color):
        # makers = 15 #标尺
        # x = 0.1 #起始坐标x
        # y = 0.1 #起始坐标y
        # x1 = 0.8#染色体长度
        # color='k'
        lw = 2  #比例线宽
        alpha = 1
        plt.axhline(y=y, xmin=x, xmax=x+x1, lw=lw, c=color, alpha=alpha)

    def d_gene(self,axes1,gff,length):
        gff = gff.sort_values(by=[2],ascending= [True])
        gff.reset_index(drop=True, inplace=True)
        print(gff)

        for index,row in gff.iterrows():
            l = abs(row[2]-row[3])
            if self.start != None:
                start = row[2] - int(self.start)
            else:
                start = row[2]

            start = (start/length)*0.8
            l = (l/length)*0.8
            # print(100*l/0.8)
            self.map2(axes1,0.1+start,0.315+(index*0.005),l,'green')


    def make_plot(self,chr_name,gff,lens):
        fig1 = plt.figure(num=1, figsize=(10, 10))  # 确保正方形在屏幕上显示一致，固定figure的长宽相等
        axes1 = fig1.add_subplot(1, 1, 1)
        x1 = 0.8
        length = lens[chr_name]['end']
        if self.start != None:
            length = length - int(self.start)
        if self.end != None:
            length = length - (lens[chr_name]['end'] - int(self.end))
        plt.xlim((0, 1))
        plt.ylim((0, 1))
        self.map1(axes1,10,0.1,0.3,x1,'#ffec47')# 长染色体0.8
        plt.text(0.45, 0.28, chr_name)
        self.d_gene(axes1,gff,length)

        plt.savefig(self.savefile,dpi=1000, bbox_inches = 'tight')

    def read_gff(self,gff,chr_list):
        gff = pd.read_csv(gff,header = None, sep='\t', comment='#')
        gff = gff.drop(gff[(gff[0] != self.chr_name)].index)
        if self.start != None:
            gff = gff.drop(gff[(gff[2] < int(self.start))].index)
        if self.end != None:
            gff = gff.drop(gff[(gff[3] > int(self.end))].index)
        return gff


    def run(self):
        # print(self.lens)
        lens,chrlist = read_lens(self.lens)
        # print(lens)
        gff = self.read_gff(self.gff,[self.chr_name])
        # print(gff)
        self.make_plot(self.chr_name,gff,lens)

