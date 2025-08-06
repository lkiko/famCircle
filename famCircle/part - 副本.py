# -*- encoding: utf-8 -*-
'''
@File        :part.py
@Time        :2021/07/11 12:44:00
@Author        :charles kiko
@Version        :1.0
@Contact        :charles_kiko@163.com
@Desc        :家族分布染色体
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
import itertools
plt.rc('font',family='Arial')



class part():
    def __init__(self, options):
        self.makers = 200
        self.a = 10
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)


    def plot_bez_inner(self, ex1x, ex1y, ex2x, ex2y):
        x = [ex1x, ex1x+((ex2x - ex1x)/3+0.1), ex2x-((ex2x - ex1x)/3+0.1), ex2x]
        y = [ex1y, ex1y+((ex2y - ex1y)/3), ex2y-((ex2y - ex1y)/3), ex2y]
        step = .01
        t = arange(0, 1+step, step)
        xt = self.Bezier(x, t)# 贝塞尔曲线
        yt = self.Bezier(y, t)
        plot(xt, yt, '-', color='r', lw=.3, alpha=0.3)#alpha 透明度

    def calculate_coef(self,p0, p1, p2, p3):
        c = 3*(p1 - p0)
        b = 3*(p2 - p1) -c
        a = p3 - p0 - c - b
        return c, b, a
    def Bezier(self,plist, t):
        p0, p1, p2, p3 = plist
        # calculates the coefficient values
        c, b, a = self.calculate_coef(p0, p1, p2, p3)
        tsquared = t**2
        tcubic = tsquared*t
        return a*tcubic + b*tsquared + c*t + p0

    def drawSpec1(self,min_,max_):#min_,max_
        # ax1 = plt.subplot(gs)
        # #波长越大ks越小
        pic = np.zeros([10,360,3])
        rgb = [getRGB(d) for d in range(400,760)]
        pic = pic+rgb
        pic0 = np.zeros([360,10,3])
        for i in range(10):
            for j in range(360):
                pic0[j][i] = pic[i][j]
        plt.imshow(pic0)
        plt.yticks(range(0,360,50),[str(format(max_ - i*((max_-min_)/8), '.4f')) for i in range(8)])      #y坐标轴ks越小颜色越红
        plt.xticks([])

    def make_plot(self,lens,gff,family,ks,family_dic,dic_motif):
        chro_lingth = lens[self.chro_name]['end']
        chro_order = lens[self.chro_name]['order']
        gene_average = chro_lingth/chro_order
        maker = 0.2
        fig1 = plt.figure(num=1, figsize=(10, 10))  # 确保正方形在屏幕上显示一致，固定figure的长宽相等
        axes1 = fig1.add_subplot(1, 1, 1)
        plt.xlim((0, 1))
        plt.ylim((0, 1))
        self.map1(axes1,maker,0.1,0.3,0.8,'k')# 长染色体0.8
        plt.text(0.45, 0.28, self.chro_name)
        plt.axis('off')
        
        y = 0.4
        for i in dic_motif.keys():
            plt.hlines(y, 0.7, 0.707,color=dic_motif[i],linewidth=2)
            plt.text(0.71, y-0.002, str(i)+"_"+str(dic_motif[i]),fontsize=4)
            y += 0.01
        # self.drawSpec1(float(self.ks_concern.split(',')[0]),float(self.ks_concern.split(',')[1]))

        self.pair_index(axes1,lens,gff,family,ks,gene_average,chro_order,chro_order,family_dic)
        plt.savefig(self.savefile,dpi=1000)

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
        for i in range(int(0.8/maker)+1):
            mx = maker*i
            plt.axvline(x=0.1+mx, ymin=y, ymax=y+(w/2.5), lw=lw, c=color, alpha=alpha)
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

    def transform_pt2(self, rad, r):
        return r*cos(rad), r*sin(rad)

    def plot_bez_Ks2(self, ex1x, ex1y, ex2x, ex2y, col, ratio):
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
        plot(xt, yt, '-', color = col, lw = 0.3)#0.1

    def cluster0(self,cluster,id1,id2):
        for i in cluster:
            if id1 in i and id2 in i:
                return True
            else:
                continue
        return False

    def clusterx(self,j,lt,gene_average,gff,dic0):
        for i in range(1,int(self.cluster)):
            if j+i in lt:
                if abs(gff[dic0[j]]['end']-gff[dic0[j+i]]['end']) < 5*gene_average:# 平均基因长度
                    return True
                else:
                    False
        return False

    def pair_index(self,axes1,lens,gff,family,ks,gene_average,chro_lingth,chro_order,family_dic):
        fam_dic = {}
        for i in family:
            if gff[i]['chr'] != self.chro_name:
                continue
            else:
                fam_dic[gff[i]['order']] = i

        # cluster = []# 基因簇
        lt = list(fam_dic.keys())
        lt.sort()# 升序
        # print(fam_dic)
        # print(lt)
        cluster = []
        cluster0 = []
        for j in lt:
            if self.clusterx(j,lt,gene_average,gff,fam_dic):# 判断目标基因相近的位置上是否存在
                cluster0.append(fam_dic[j])
            else:
                cluster0.append(fam_dic[j])
                cluster.append(cluster0)
                cluster0 = []
        # print(cluster)
        clusterx = {}
        for i in cluster:
            for j in range(len(i)):
                id0 = i[j]
                clusterx[id0] = j
                x1 = (0.8/chro_lingth)*int(self.cluster)
                print(x1)
                self.map2(axes1, ((gff[id0]['order']/chro_lingth)*0.8)+0.1, (j)*0.005+0.32+(j)*0.001, x1, family_dic[id0])
        ratio = 0.5
        put1 = []
        for lt in cluster:
            if len(lt) <= 1:
                continue
            pair = list(itertools.permutations(lt,2))
            for i in pair:
                id1,id2 = i[0],i[1]
                if id1 in family and id2 in family and self.cluster0(cluster,id1,id2):
                    if self.parameter == 'True':
                        if str(id1)+'_'+str(id2) in ks.keys():
                            ks_v = ks[str(id1)+'_'+str(id2)]
                        elif str(id2)+'_'+str(id1) in ks.keys():
                            ks_v = ks[str(id2)+'_'+str(id1)]
                        else:
                            continue
                        if ks_v >= float(self.ks_concern.split(',')[0]) and ks_v <= float(self.ks_concern.split(',')[1]):
                            # col = return_col(ks_v,float(self.ks_concern.split(',')[0]),float(self.ks_concern.split(',')[1]))
                            put1.append([id1,id2])
                            continue
                        else:
                            col = "#d4dcda"
                        self.plot_bez_Ks2(((gff[id1]['order']/chro_lingth)*0.8)+0.101, (clusterx[id1]) * 0.006+0.32, ((gff[id2]['order']/chro_lingth)*0.8)+0.101, (clusterx[id2])*0.006+0.32, col, ratio)
                    else:
                        self.plot_bez_Ks2(((gff[id1]['order']/chro_lingth)*0.8)+0.101, (clusterx[id1]) * 0.006+0.32, ((gff[id2]['order']/chro_lingth)*0.8)+0.101, (clusterx[id2])*0.006+0.32, 'red', ratio)

        if self.parameter == 'True':
            for i in put1:
                id1,id2 = i[0],i[1]
                if id1 in family and id2 in family and self.cluster0(cluster,id1,id2):
                    if str(id1)+'_'+str(id2) in ks.keys():
                        ks_v = ks[str(id1)+'_'+str(id2)]
                    elif str(id2)+'_'+str(id1) in ks.keys():
                        ks_v = ks[str(id2)+'_'+str(id1)]
                    else:
                        continue
                    col = return_col(ks_v,float(self.ks_concern.split(',')[0]),float(self.ks_concern.split(',')[1]))
                    self.plot_bez_Ks2(((gff[id1]['order']/chro_lingth)*0.8)+0.101, (clusterx[id1]) * 0.006+0.32, ((gff[id2]['order']/chro_lingth)*0.8)+0.101, (clusterx[id2])*0.006+0.32, col, ratio)

    def run(self):
        lens,chr_list = read_lens(self.lens)# 染色体字典
        gff = read_gff(self.gff)# 基因字典
        family,family_dic,dic_motif = read_family1(self.genefamily)# 家族列表
        for i in family:# 去除不在染色体上的家族成员
            if i not in gff.keys():
                family.remove(i)
        ks = read_ks0(self.ks)# ks字典
        self.make_plot(lens,gff,family,ks,family_dic,dic_motif)
