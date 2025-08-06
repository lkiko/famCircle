# -*- encoding: utf-8 -*-
'''
@File        :outer.py
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
plt.rc('font',family='Arial')


class part_out():
    def __init__(self, options):
    ##### circle parameters
        self.block = 0
        self.dpi = 600
        self.start = 0
        self.end = -1
        self.gap_ratio = 20 #gaps between chromosome circle, chr:gap = 4: 1
        self.radius = 0.3
        self.chrscize = 0.01
        self.block_scize= 0.004 #block scize
        self.blockhight = 0.005
        self.blockthick = 0.007 #0.006
        self.distance_down = 100000000
        self.distance_up = 10000000000
        self.ks_file_type = 'ks'
        self.model = 'gradual'
        self.node = '0,0.1,0.2,0.5,1'
        self.distance_boolean = "False"
        self.start_list = {}
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)

    def ksrun(self):
        lens,chr_list = read_lens(self.lens)# 染色体字典
        gff = read_gff(self.gff,chr_list)# 基因字典
        family,family_dic,dic_motif = read_family1(self.genefamily)# 家族列表
        family = drop_none(family,gff)
        ks,ks_l = read_ks1(self.ks,self.ks_file_type,family)# ks字典
        return lens,chr_list,gff,ks,family,family_dic,dic_motif,ks_l

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

    def plot_arc_block(self, start, radius,col0):# block
        t = arange(start, start+self.block_scize, pi/720.)
        x,y = radius * cos(t), radius*sin(t)
        x1, y1 = (radius-self.blockthick) * cos(t), (radius-self.blockthick) * sin(t)
        plot(x, y, col0, linewidth=3, alpha=0.9)

    def zj(self,lens,chr_list):
        fullchrolen = sum([lens[i]['end'] for i in lens.keys()])
        fullgene = sum([lens[i]['order'] for i in lens.keys()])
        gene_average = int(fullchrolen/fullgene)
        # print("gene_average     ",gene_average)
        chr_number = len(lens.keys()) # total number of chromosomes
        GAP = fullchrolen/float(self.gap_ratio)/chr_number # gap size in base pair
        total_size = fullchrolen + chr_number * GAP # base pairs 
        # print('total_size',total_size)
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
        plot(xt, yt, '-', color=cl, lw=lw0, alpha=0.5)#alpha 透明度

    def plot_bez_Ks2(self, rad1, r1, rad2, r2, col, ratio):
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

    def run(self):
        self.radius_a =  float(self.radius)
        self.radius_b = self.radius_a + float(self.chrscize) #.33, .335   # 半径控制参数 
        self.sm_radius=(self.radius_b-self.radius_a)/2 #telomere capping
        if (os.path.exists(self.savecsv)):
            os.remove(self.savecsv)
        lens,chr_list,gff,ks,family,family_dic,dic_motif,ks_l = self.ksrun()
        chr_list = [self.chr_name]
        lens0={}
        lens0[self.chr_name] = lens[self.chr_name]
        del lens
        lens = lens0
        total_size,gene_average = self.zj(lens,chr_list)
        stop_list = {}

        fig = plt.figure(figsize=(10, 11))
        gs = gridspec.GridSpec(30, 30)#, width_ratios=[10, 1]
        ax1 = fig.add_subplot(gs[1:25,1:25], aspect=1)

        for i in self.start_list.keys():
            stop_list[i] = self.start_list[i] + lens[i]['end']
        for i in chr_list:
            start,stop = self.start_list[i],stop_list[i]
            start, stop = self.to_radian(start, total_size), self.to_radian(stop, total_size)
            # shaft
            x,y = plot_arc(start, stop,self.radius_a,self.chrscize)
            plt.fill(x, y, facecolor='red',alpha=.7)

        # 绘制关系，并查询outer
        chro_order = lens[self.chr_name]['order']
        if self.end == -1:
            end = chro_order
            chro_order = chro_order - int(self.start)
        else:
            end = int(self.end)
            chro_order = int(self.end) - int(self.start)

        dicx0 = {}
        fam_list = []
        for i in family:
            if gff[i]['chr'] != self.chr_name or gff[i]['order'] < int(self.start) or gff[i]['order'] > int(end):# 判断家族成员是否在绘制区
                continue
            dicx0[gff[i]['order']] = i# 创建相对位置-> 基因
            fam_list.append(i)

        cluster = []# 基因簇
        lt = list(dicx0.keys())
        lt.sort()# 相对位置升序

        cluster0 = []
        for i in range(len(lt)):# j为相对位置
            j = lt[i]
            if clusterx_pro(i,j,lt,gene_average,gff,dicx0,self.cluster,self.magnification):# 相对位置中存在串联单位
                cluster0.append(dicx0[j])# 添加到子基因簇
            else:
                if len(cluster0) > 0:
                    cluster.append(cluster0)
                cluster0 = []
                cluster0.append(dicx0[j])# 子基因簇归位
        cluster.append(cluster0)
        # print(cluster0)

        lt = ['number','chr', 'gene_list']
        with open(self.savecsv, "a", newline='', encoding='utf-8') as file:
            writer = csv.writer(file ,delimiter=',')
            writer.writerow(lt)
        fam_clu = {}
        num_cl = 0
        for i in cluster:
            num_cl += 1
            lt = [num_cl,gff[i[0]]['chr']] + [','.join(i)]
            with open(self.savecsv, "a", newline='', encoding='utf-8') as file:
                writer = csv.writer(file ,delimiter=',')
                writer.writerow(lt)
            # 输出簇编号
            if self.parameter == 'True':
                posi = gff[i[0]]['end'] + self.start_list[gff[i[0]]['chr']]
                start = self.to_radian(posi, total_size)
                label_x, label_y = self.rad_to_coord(start, self.radius_b*0.965)#1.2
                #print label_x, label_y
                text(label_x, label_y, str(num_cl), horizontalalignment="center", verticalalignment="center", fontsize = 3.5, color = 'black')
            for j in range(len(i)):
                fam_clu[i[j]] = j
                posi = gff[i[j]]['end'] + self.start_list[gff[i[j]]['chr']]
                # start = self.to_radian(posi, total_size)
                start, stop = self.to_radian(posi, total_size),self.to_radian(posi+(total_size*self.block_scize), total_size)
                col0 = family_dic[i[j]]
                # self.plot_arc_block(start, self.radius_b + (j+1) * self.blockthick,col0)# 基因
                # col0 = family_dic[i[j]]
                xb,yb = plot_block_s(start, stop ,self.blockhight,self.radius_b + (j+1) * self.blockthick)
                plt.fill(xb, yb, facecolor=col0,alpha=1)


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
        del ks,ks_m,df_clear
        lt = df[0].to_list()
        ks = df.set_index(0)[1].to_dict()

        print("分类ks")
        inner = []
        outer = []
        for j in trange(len(lt)):
            i = lt[j]
            id1,id2 = i.split("%")[0],i.split("%")[1]
            if id1 in fam_list and id2 in fam_list and self.cluster0(cluster,id1,id2):
                outer.append([id1,id2])# 簇内
            elif id1 in fam_list and id2 in fam_list:
                inner.append([id1,id2])# 簇间

        output = open("part_out_"+str(self.chr_name)+"_"+str(self.ks_concern.split(',')[0])+"_"+str(self.ks_concern.split(',')[1])+".outer","w")
        output.write("id1\tid2\tchro1\tchro2\tstart1\tstart2\tks_v\tdistance\tTYPE\n")

        y = -0.35
        for i in dic_motif.keys():
            plt.hlines(y, 0.3, 0.307,color=dic_motif[i],linewidth=2)
            plt.text(0.31, y-0.002, str(i)+"_"+str(dic_motif[i]),fontsize=6)
            y += 0.01

        ratio = 0.5

        print("绘制基因簇内部 ",len(outer))
        for j in trange(len(outer)):
            i = outer[j]
            ks_v = ks['%'.join(i)]
            id1,id2 = i[0],i[1]
            if ks_v >= float(self.ks_concern.split(',')[0]) and ks_v <= float(self.ks_concern.split(',')[1]):
                col = return_col1(self.model,self.node,ks_v,float(self.ks_concern.split(',')[0]),float(self.ks_concern.split(',')[1]))
                posi1 = gff[id1]['end'] + self.start_list[gff[id1]['chr']]+ (total_size*self.block_scize/2)
                start1 = self.to_radian(posi1, total_size)
                posi2 = gff[id2]['end'] + self.start_list[gff[id2]['chr']]+ (total_size*self.block_scize/2)
                start2 = self.to_radian(posi2, total_size)
                output.write(str(id1)+"\t"+str(id2)+"\t"+str(gff[id1]['chr'])+"\t"+str(gff[id2]['chr'])+"\t"+\
                    str(gff[id1]['end'])+"\t"+str(gff[id2]['end'])+"\t"+str(ks_v)+"\t"+str(abs(gff[id2]['end']-gff[id1]['end']))+"\tgene_cluster\n")
                self.plot_bez_Ks2(start1, (self.radius_b + (fam_clu[id1]+1) * self.blockthick+(self.blockhight/2)), start2, self.radius_b + (fam_clu[id2]+1) * self.blockthick+(self.blockhight/2), col, ratio)

        print("绘制基因簇之间 ",len(inner))
        dic01 = {}
        for j in trange(len(inner)):
            i = inner[j]
            # print(i)
            ks_v = ks['%'.join(i)]
            if ks_v >= float(self.ks_concern.split(',')[0]) and ks_v <= float(self.ks_concern.split(',')[1]):
                # print(i)
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
            pos1,pos2 = gff[id1]['end'],gff[id2]['end']
            order1,order2 = gff[id1]['order'],gff[id2]['order']
            chro1,chro2 = gff[id1]['chr'],gff[id2]['chr']
            if self.distance_boolean == "True":
                if abs(pos1-pos2) < self.distance_down or abs(pos1-pos2) > self.distance_up:
                    continue
            ks_v = dic01[str(id1) + '%' + str(id2)]
            col = return_col1(self.model,self.node,ks_v,float(self.ks_concern.split(',')[0]),float(self.ks_concern.split(',')[1]))
            output.write(str(id1)+"\t"+str(id2)+"\t"+str(chro1)+"\t"+str(chro2)+"\t"+str(pos1)+"\t"+str(pos2)\
                +"\t"+str(ks_v)+"\t"+str(abs(pos1-pos2))+"\tbetween_gene_cluster\n")
            self.plot_bez_inner((chro1, pos1, self.radius_a*0.95), (chro2, pos2, self.radius_a*0.95), col, total_size,1)#绘制内部共线性

        output.close()

        plt.axis('off')
        if self.model == 'gradual':
            drawSpec(float(self.ks_concern.split(',')[0]),float(self.ks_concern.split(',')[1]),gs[:,25:])
        elif self.model == 'jump':
            drawSpec_1(self.node,gs[:,25:])
        savefig(self.savefile, dpi = int(self.dpi), bbox_inches = 'tight')
        sys.exit(0)