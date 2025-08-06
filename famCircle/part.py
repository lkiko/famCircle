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
import gc
from matplotlib import gridspec
from tqdm import trange
plt.rc('font',family='Arial')

class part():
    def __init__(self, options):
        self.makers = 200
        self.start = float ("-inf")
        self.end = float ("inf")
        self.block_scize= 0.004 #block scize
        self.blockhight = 0.005
        self.inter_cluster = "True"
        self.ks_file_type = 'ks'
        self.distance_down = 100000000
        self.distance_up = 10000000000
        self.distance_boolean = "False"
        self.model = 'gradual'
        self.node = '0,0.1,0.2,0.5,1'
        self.parameter = 'True'
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

    def make_plot(self,lens,gff,family,ks,family_dic,dic_motif,ks_l):
        chro_lingth = lens[self.chr_name]['end']
        chro_order = lens[self.chr_name]['order']
        gene_average = chro_lingth/chro_order
        # print('gene_average is ',gene_average)
        maker = 0.2
        fig = plt.figure(figsize=(10, 10))
        gs = gridspec.GridSpec(30, 30)#, width_ratios=[10, 1]
        axes1 = fig.add_subplot(gs[1:25,:])#, aspect=1

        plt.xlim((0, 1))
        plt.ylim((0, 1))

        if self.inter_cluster == "True":# 簇间关系开关
            start_d = 0.48
        else:
            start_d = 0.05
        self.map1(axes1,maker,0.05,start_d,0.9,'k')# 染色体0.9

        plt.text(0.45, 0, self.chr_name)
        plt.axis('off')

        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())

        y = 0.7
        for i in dic_motif.keys():# 标签绘制
            plt.hlines(y, 0.7, 0.707,color=dic_motif[i],linewidth=2)
            plt.text(0.71, y-0.002, str(i)+"_"+str(dic_motif[i]),fontsize=6)
            y += 0.015

        self.pair_index(axes1,lens,gff,family,ks,gene_average,chro_lingth,chro_order,family_dic,ks_l)# 主要绘制
        if self.model == 'gradual':
            drawSpec1(float(self.ks_concern.split(',')[0]),float(self.ks_concern.split(',')[1]),gs[25:,:])
        elif self.model == 'jump':
            drawSpec11(self.node,gs[25:,:])
        plt.savefig(self.savefile,dpi=1500, bbox_inches = 'tight')

    def map1(self,axes1,maker,x,y,x1,color):
        # makers = 15 #标尺
        # x = 0.1 #起始坐标x
        # y = 0.1 #起始坐标y
        # x1 = 0.9#染色体长度
        # color='k'
        lw = 1  #比例线宽
        y1 = 0.001
        w = 0.01
        alpha = 1
        plt.axhline(y=y, xmin=x, xmax=x+x1, lw=lw, c=color, alpha=alpha)
        plt.axhline(y=y+w+y1, xmin=x, xmax=x+x1, lw=lw, c=color, alpha=alpha)

        for i in range(5+1):
            mx = (0.9/5)*i
            plt.axvline(x=0.05+mx, ymin=y, ymax=y+(w/2.5), lw=lw, c=color, alpha=alpha)
        if self.end != None:
            o = self.end
            xt = 0.02+mx
            plt.text(xt, y-0.015, str(o),fontsize=6)
            return 0
        else:
            o = '+'
            xt = 0.045+mx
            plt.text(xt, y-0.015, str(o),fontsize=6)

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
        # x1 = 0.9#染色体长度
        # color='k'
        lw = 2  #比例线宽
        alpha = 0.8
        # y = y - 0.00223
        # plt.axvline(x=x, ymin= y ,ymax = y + 0.005,lw = 0.3, c=color, alpha=alpha)
        # print(x,x1,x+x1)
        plt.axhline(y=y, xmin=x, xmax=x + x1, lw=lw, c=color, alpha=alpha)

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

    def clusterx(self,index,j,lt,gene_average,gff,dic0):
        for i in range(1,int(self.cluster)+1):
            if j+i in lt or j-i in lt:
                dow = j-i
                up = j+i
                if dow in lt and abs(gff[dic0[j]]['end']-gff[dic0[j-i]]['end']) < int(self.cluster)*gene_average*int(self.magnification):
                    # print(gff[dic0[j]]['end'],gff[dic0[j-i]]['end'],abs(gff[dic0[j]]['end']-gff[dic0[j-i]]['end']),int(self.cluster)*gene_average*int(self.magnification))
                    return True
                elif up in lt and abs(gff[dic0[j]]['end']-gff[dic0[j+i]]['end']) < int(self.cluster)*gene_average*int(self.magnification):# 平均基因长度
                    if index != 0 and abs(gff[dic0[j]]['end']-gff[dic0[lt[index-1]]]['end']) < int(self.cluster)*gene_average*int(self.magnification):
                        # print(gff[dic0[j]]['end'],gff[dic0[lt[index-1]]]['end'],abs(gff[dic0[j]]['end']-gff[dic0[lt[index-1]]]['end']),int(self.cluster)*gene_average*int(self.magnification))
                        return True
                    else: 
                        return False
                else:
                    False
        return False

    def get_gff(self,gff,name):
        return gff[name]['chr']

    def pair_index(self,axes1,lens,gff,family,ks,gene_average,chro_lingth,chro_order,family_dic,ks_l):
        if self.inter_cluster == "True":
            start_d = 0.5
        else:
            start_d = 0.07
        start = int(max(float(self.start),0))
        end =  int(min(float(self.end),chro_lingth))
        chro_lingth = abs(end - start)
        # if self.end == None:
        #     end = chro_lingth
        #     chro_lingth = chro_lingth - int(self.start)
        # else:
        #     end = int(self.end)
        #     chro_lingth = int(self.end) - int(self.start)

        fam_dic = {}
        fam_list = []
        for i in family:
            if gff[i]['chr'] != self.chr_name or gff[i]['end'] < int(start) or gff[i]['end'] > int(end):
                continue
            else:
                fam_dic[gff[i]['order']] = i
                fam_list.append(i)

        lt = list(fam_dic.keys())
        lt.sort()# 升序
        # print(lt)
        cluster = []
        cluster0 = []
        for i in range(len(lt)):
            j = lt[i]
            if clusterx_pro(i,j,lt,gene_average,gff,fam_dic,self.cluster,self.magnification):# 判断目标基因相近的位置上是否存在
                cluster0.append(fam_dic[j])
            else:
                if len(cluster0) > 0:
                    cluster.append(cluster0)
                cluster0 = []
                cluster0.append(fam_dic[j])# 子基因簇归位
        cluster.append(cluster0)

        width = (0.9/(chro_lingth/gene_average))#*(int(self.cluster))
        clusterx = {}

        lt = ['number','chr', 'gene_list']
        with open(self.savecsv, "a", newline='', encoding='utf-8') as file:
            writer = csv.writer(file ,delimiter=',')
            writer.writerow(lt)
        num_cl = 0
        for i in cluster:
            for j in range(len(i)):
                id0 = i[j]
                clusterx[id0] = j
                # posi = gff[i[0]]['end']-int(start)
                posi = gff[id0]['end']-int(start)
                print(posi)
                # continue
                start, stop = (posi/chro_lingth)*0.9+0.05,((posi+(chro_lingth*self.block_scize))/chro_lingth)*0.9+0.05
                xb,yb = plot_block_s_line(start, stop ,self.blockhight,(j)*0.005+start_d+(j)*0.001)
                plt.fill(xb, yb, facecolor=family_dic[id0],alpha=1)
                # self.map2(axes1, (((gff[id0]['end']-int(start))/chro_lingth)*0.9)+0.05, \
                #     (j)*0.005+start_d+(j)*0.001, width, family_dic[id0])

            num_cl += 1
            lt = [num_cl,gff[i[0]]['chr']] + [','.join(i)]
            with open(self.savecsv, "a", newline='', encoding='utf-8') as file:
                writer = csv.writer(file ,delimiter=',')
                writer.writerow(lt)

            if self.parameter == 'True':
                posi = (((gff[i[0]]['end']-int(start))/chro_lingth)*0.9)+0.05
                text(posi, start_d-0.026, str(num_cl), horizontalalignment="center", verticalalignment="center", fontsize = 3.5, color = 'black')

        ratio = 0.5
        put1 = []
        ks_m = pd.DataFrame(list(ks.items()))

        ks_m[2] = ks_m[0].apply(qiepian0)
        ks_m[3] = ks_m[0].apply(qiepian1)
        ks_m[4] = ks_m[2].apply(lambda n: gff[n]['chr'])
        ks_m[5] = ks_m[3].apply(lambda n: gff[n]['chr'])

        df_clear = ks_m.drop(ks_m[(ks_m[4] != self.chr_name) | (ks_m[5] != self.chr_name)].index) 
        df = df_clear[[0,1]]# 保留两列
        if self.ks_file_type == 'ks':
            df = df.sort_values(by=1, ascending=False)#降序
        elif self.ks_file_type == 'blast':
            df = df.sort_values(by=1, ascending=True)#升序
        lt = df[0].to_list()
        del ks,ks_m,df_clear
        ks = df.set_index(0)[1].to_dict()
        # print(ks)
        print("分类ks")
        inner = []# 簇间
        outer = []# 簇内部
        for j in trange(len(lt)):
            i = lt[j]
            id1,id2 = i.split("%")[0],i.split("%")[1]
            if id1 in fam_list and id2 in fam_list and self.cluster0(cluster,id1,id2):
                outer.append([id1,id2])# 簇内
            elif id1 in fam_list and id2 in fam_list:
                inner.append([id1,id2])# 簇间

        print("基因簇内部    ",len(outer))
        output = open("part_"+str(self.chr_name)+"_"+str(self.ks_concern.split(',')[0])+"_"+str(self.ks_concern.split(',')[1])+".outer","w")
        output.write("id1\tid2\tchro1\tchro2\tstart1\tstart2\tks_v\tdistance\tTYPE\n")
        for j in trange(len(outer)):
            i = outer[j]
            ks_v = ks['%'.join(i)]
            id1,id2 = i[0],i[1]
            chro1,chro2 = gff[id1]['chr'],gff[id2]['chr']
            pos1,pos2 = gff[id1]['end'],gff[id2]['end']
            # if self.parameter == 'True':
            if ks_v >= float(self.ks_concern.split(',')[0]) and ks_v <= float(self.ks_concern.split(',')[1]):
                output.write(str(id1)+"\t"+str(id2)+"\t"+str(chro1)+"\t"+str(chro2)+"\t"+str(pos1)+"\t"+str(pos2)\
                    +"\t"+str(ks_v)+"\t"+str(abs(pos1-pos2))+"\tgene_cluster\n")
                col = return_col1(self.model,self.node,ks_v,float(self.ks_concern.split(',')[0]),float(self.ks_concern.split(',')[1]))
                # self.plot_bez_Ks2((((gff[id1]['end']-int(start))/chro_lingth)*0.9)+0.050, \
                #     (clusterx[id1]) * 0.006+start_d, (((gff[id2]['end']-int(start))/chro_lingth)*0.9)+0.050, \
                #     (clusterx[id2])*0.006+start_d, col, ratio)
                self.plot_bez_Ks2((((gff[id1]['end']-int(start)+(chro_lingth*self.block_scize/2))/chro_lingth)*0.9+0.05), \
                    (clusterx[id1]) * 0.006+start_d+(self.blockhight/2), (((gff[id2]['end']-int(start)+(chro_lingth*self.block_scize/2))/chro_lingth)*0.9+0.05), \
                    (clusterx[id2])*0.006+start_d+(self.blockhight/2), col, ratio)
        if self.inter_cluster == "True":
            print("基因簇之间    ",len(inner))
            dic01 = {}
            for j in trange(len(inner)):
                i = inner[j]
                ks_v = ks['%'.join(i)]
                if ks_v >= float(self.ks_concern.split(',')[0]) and ks_v <= float(self.ks_concern.split(',')[1]):
                    dic01['%'.join(i)] = ks_v
            L = list(dic01.items())
            if float(self.ks_concern.split(',')[1]) < 10:
                L.sort(key=lambda x:x[1],reverse=True)# 从大到小适合ks
            else:
                L.sort(key=lambda x:x[1],reverse=False)# 从小到大适合blast
            lt01 = []
            for i in L:
                lt01.append([i[0].split("%")[0],i[0].split("%")[1]])
            for j in trange(len(lt01)):
                i = lt01[j]
                id1,id2 = i[0],i[1]
                ks_v = dic01[str(id1) + '%' + str(id2)]
                chro1,chro2 = gff[id1]['chr'],gff[id2]['chr']
                pos1,pos2 = gff[id1]['end'],gff[id2]['end']
                if self.distance_boolean == "True":
                    if abs(pos1-pos2) < self.distance_down or abs(pos1-pos2) > self.distance_up:
                        continue
                output.write(str(id1)+"\t"+str(id2)+"\t"+str(chro1)+"\t"+str(chro2)+"\t"+str(pos1)+"\t"+str(pos2)\
                    +"\t"+str(ks_v)+"\t"+str(abs(pos1-pos2))+"\tbetween_gene_cluster\n")
                col = return_col1(self.model,self.node,ks_v,float(self.ks_concern.split(',')[0]),float(self.ks_concern.split(',')[1]))
                self.plot_bez_Ks2((((gff[id1]['end']-int(start))/chro_lingth)*0.9)+0.050+(width/2), \
                    start_d-0.035, (((gff[id2]['end']-int(start))/chro_lingth)*0.9)+0.050+(width/2), \
                    start_d-0.035, col, ratio)

    def run(self):
        lens,chr_list = read_lens(self.lens)# 染色体字典
        gff = read_gff(self.gff,chr_list)# 基因字典
        family,family_dic,dic_motif = read_family1(self.genefamily)# 家族列表
        family = drop_none(family,gff)
        ks,ks_l = read_ks1(self.ks,self.ks_file_type,family)# ks字典
        self.make_plot(lens,gff,family,ks,family_dic,dic_motif,ks_l)
