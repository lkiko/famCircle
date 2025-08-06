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

class sline():
    def __init__(self, options):
        self.makers = 200
        self.alpha_f = .7
        self.linewidth = 0.2
        self.class1 = False
        self.dpi = 1000
        self.plot_h = 10
        self.plot_w = 10
        self.block = 6
        self.block_gep = 200
        self.position = 'order'
        self.fucus_lw = 2
        self.type = 'None'
        self.gap = 4
        self.h1,self.h2,self.h = 0.1,.9,.02
        self.blast_reverse = 'False'
        self.chr_name_size = '12'
        self.focus_list = 'None'
        self.focus_p = 'True'
        self.focusc = 'red'
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)
        self.block = int(self.block)
        self.block_gep = int(self.block_gep)
        self.gap = float(self.gap)
        self.fucus_lw = float(self.fucus_lw)
        self.h1 = float(self.h1)
        self.h2 = float(self.h2)
        self.linewidth = float(self.linewidth)
        self.alpha_f = float(self.alpha_f)
        self.h = float(self.h)
        self.plot_h = int(self.plot_h)
        self.plot_w = int(self.plot_w)

        self.colord = {}
        if ';' in self.cols:
            lt = self.cols.split(';')
            for i in range(len(lt)):
                if ':' in lt[i]:
                    lt1 = lt[i].split(',')
                    for j in range(len(lt1)):
                        self.colord[lt1[j].split(':')[0]] = lt1[j].split(':')[1]
                else:
                    self.colord[i] = lt[i]
        else:
            lt = self.cols.split(',')
            for j in range(len(lt)):
                self.colord[j] = lt[j]
        print(self.colord)

    def drow_coli0(self,s1hpop,s2hpop,h1,h2,p):
        h = self.h
        if p:
            xt,yt = plot_bezier_curve(s1hpop,h1-(h)*0.1,s2hpop,h2+(h)*1.1)
        else:
            xt,yt = plot_bezier_curve(s1hpop,h1+(h)*1.1,s2hpop,h2-(h)*0.1)

        return xt,yt

    def drow_coli(self,s1hpop,s2hpop,s1dpop,s2dpop,h1,h2,p):
        h = self.h
        x,y = [],[]
        if p:
            x.append(s1hpop)
            y.append(h1-(h)*0.1)
            x.append(s1dpop)
            y.append(h1-(h)*0.1)
            xt,yt = plot_bezier_curve(s1dpop,h1-(h)*0.1,s2dpop,h2+(h)*1.1)
        else:
            x.append(s1hpop)
            y.append(h1+(h)*1.1)
            x.append(s1dpop)
            y.append(h1+(h)*1.1)
            xt,yt = plot_bezier_curve(s1dpop,h1+(h)*1.1,s2dpop,h2-(h)*0.1)
        x = x + [i for i in xt]
        y = y + [i for i in yt]
        if p:
            x.append(s2dpop)
            y.append(h2+(h)*1.1)
            x.append(s2hpop)
            y.append(h2+(h)*1.1)
            xt,yt = plot_bezier_curve(s2hpop,h2+(h)*1.1,s1hpop,h1-(h)*0.1)
        else:
            x.append(s2dpop)
            y.append(h2-(h)*0.1)
            x.append(s2hpop)
            y.append(h2-(h)*0.1)
            xt,yt = plot_bezier_curve(s2hpop,h2-(h)*0.1,s1hpop,h1+(h)*1.1)
        x = x + [i for i in xt]
        y = y + [i for i in yt]
        return x,y

    def run(self):
        lenss,chrlists = [],[]
        for file in self.lens.split(','):
            lens1,chrlist1 = read_lens(file)
            # print(lens1,chrlist1)
            lenss.append(lens1)
            chrlists.append(chrlist1)
        # exit()
        chr_color = {}
        if self.chr_name != 'all':
            chrlists = []
            for s in self.chr_name.split(';'):
                if ':' in s:
                    lt0 = []
                    list0 = s.split(',')
                    for i in list0:
                        if ':' not in i:
                            lt0.append(i)
                        else:
                            lt = i.split(':')
                            chr_color[lt[0]] = lt[1]
                            lt0.append(lt[0])
                    chrlists.append(lt0)
                else:
                    chrlists.append(s.split(','))

        gff = None
        for i in range(len(lenss)):
            print(i)
            print(chrlists[i])
            gfff = read_gff(self.gffs.split(',')[i],chrlists[i])
            # print(gfff)
            if gff == None:
                gff = gfff
            else:
                gff = dict( gff, **gfff)

        if self.focus_list != 'None':
            focus_list_pd = pd.read_csv(self.focus_list,header = None, sep='\t', comment='#')
            if len(focus_list_pd.columns) == 1:
                focus_list = list(focus_list_pd[0].to_list())
                focus_pair = []
                focus_list_pd[1] = self.focusc
                focus_listc = dict(zip(focus_list_pd[0],focus_list_pd[1]))
            elif len(focus_list_pd.columns) == 3:
                focus_list = list(set(list(focus_list_pd[0].to_list()) + list(focus_list_pd[1].to_list())))
                focus_pair = [(a,b) for a,b in zip(focus_list_pd[0].to_list(),focus_list_pd[1].to_list())]
                # focus_pair_mapping = [a+'lanmeifang'+b for a,b in zip(focus_list_pd[0].to_list(),focus_list_pd[1].to_list())]
                focus_list_pd[3] = focus_list_pd[0] + 'lanmeifang' + focus_list_pd[1]
                # focus_list_pd[1] = self.focusc
                focus_listc = dict(zip(focus_list_pd[3],focus_list_pd[2]))
            elif len(focus_list_pd.columns) == 2 and focus_list_pd.iloc[0, 1] in gff.keys():
                focus_list = list(set(list(focus_list_pd[0].to_list()) + list(focus_list_pd[1].to_list())))
                focus_pair = [(a,b) for a,b in zip(focus_list_pd[0].to_list(),focus_list_pd[1].to_list())]
                focus_list_pd[2] = focus_list_pd[0] + 'lanmeifang' + focus_list_pd[1]
                focus_list_pd[3] = self.focusc
                focus_listc = dict(zip(focus_list_pd[2],focus_list_pd[3]))
            elif len(focus_list_pd.columns) == 2 and focus_list_pd.iloc[0, 1] not in gff.keys():
                focus_list = list(focus_list_pd[0].to_list())
                focus_pair = []
                focus_listc = dict(zip(focus_list_pd[0],focus_list_pd[1]))
            else:
                print('The focus_list file format supports the following structures: gene1\ngene2\ngene3\n or \ngene1\tgene2\ngene3\tgene4\n or \ngene1\tcolor\ngene2\tcolor\n or \ngene1\tgene2\tcolor\ngene3\tgene4\tcolor.')

        # 只保留需要的染色体共线性
        lt = self.genepairs.split(',')
        colineartlys,focus_list2 = [],[]
        for i in range(len(lt)):
            if self.type != 'line':
                if self.focus_list == 'None':
                    colineartly = readcolineartly(lt[i],self.block,self.block_gep,self.genepairsfile_type,gff,chrlists[i]+chrlists[i+1])
                    colineartlys.append(colineartly)
                else:
                    colineartly,focus_list20 = readcolineartly1(lt[i],self.block,self.block_gep,self.genepairsfile_type,gff,chrlists[i]+chrlists[i+1],focus_list,self.focus_p)
                    colineartlys.append(colineartly)
                    focus_list2.append(focus_list20)
            else:
                # print(chrlists[i]+chrlists[i+1])
                colineartlys.append(read_coliearity_sline(lt[i],self.block,self.genepairsfile_type,chrlists[i]+chrlists[i+1]))
        # print(colineartlys)
        sumchrs = []
        for i in range(len(lenss)):
            sumchr = sum([lenss[i][j][self.position] for j in chrlists[i]])
            sumchrs.append(sumchr)

        # total_size = ((self.gap+1)*max(sumchrs)/self.gap)# 160508864
        total_size = ((self.gap+1)*160508864/self.gap)# 160508864
        step = 0.99/total_size
        ns = []
        for lt in chrlists:
            ns.append(len(lt))
        se_dics = []
        for i in range(len(ns)):
            se_dic = {}
            GAP1 = (total_size-sumchrs[i])/(ns[i]+1)
            start = GAP1
            for j in range(ns[i]):
                i0 = chrlists[i][j]
                dic = {}
                dic['start'] = start
                dic['end'] = start + lenss[i][i0][self.position]
                start = start + lenss[i][i0][self.position] + GAP1
                se_dic[i0] = dic
            se_dics.append(se_dic)

        fig, ax = plt.subplots(figsize=(self.plot_h, self.plot_w))
        path = get_path()
        font_path = os.path.join(path, 'example/arial.ttf')
        from matplotlib.font_manager import FontProperties
        font_prop = FontProperties(fname=font_path)
        ax.set_aspect('equal')

        h0 = self.h2#起始高度
        mh = (self.h2-self.h1)/len(ns)
        for i in range(len(ns)):
            for chro in chrlists[i]:
                # print(chro)
                xmin,xmax = se_dics[i][chro]['start']*step,se_dics[i][chro]['end']*step
                x,y = line_chr(h0,xmin,xmax,self.h)
                plt.text((xmin+xmax)/2, h0-0.02,chro,ha='center',va='center', fontproperties=font_prop,fontsize=int(self.chr_name_size))
                if chro in self.colord.keys():
                    color = self.colord[chro]
                else:
                    color = self.colord[i]
                plt.fill(x,y, facecolor=color,alpha=.7)# self.cols.split(',')[i]
            h0 -= mh

        h0 = self.h2#起始高度
        for num in range(len(ns)-1):
            if self.type != 'line':
                for s1h,s2h,s1d,s2d in colineartlys[num]:
                    if self.blast_reverse == "True":
                        s1h,s2h = s2h,s1h
                        s1d,s2d = s2d,s1d

                    s1hpop,s2hpop,s1dpop,s2dpop = gff[s1h][self.position],gff[s2h][self.position],gff[s1d][self.position],gff[s2d][self.position]
                    s1hchr,s2hchr,s1dchr,s2dchr = gff[s1h]['chr'],gff[s2h]['chr'],gff[s1d]['chr'],gff[s2d]['chr']
                    if s1hchr in chr_color.keys():
                        chro_color = chr_color[s1hchr]
                    elif s2hchr in chr_color.keys():
                        chro_color = chr_color[s2hchr]
                    else:
                        chro_color = self.linecs.split(',')[num]


                    if s1hchr not in (chrlists[num]+chrlists[num+1]) or s2hchr not in (chrlists[num]+chrlists[num+1]):# 排除不在绘制范围内的共线性
                        continue

                    if s1hchr != s2hchr:# 如果当前block在不同染色体上
                        # print(s1hchr,s2hchr)
                        if self.genepairsfile_type == 'MCScanX':
                            if s1hchr in chrlists[num] and s2dchr in chrlists[num+1]:# 共线性第一列就是第二排
                                # s1h,s2h = s2h,s1h
                                # s1d,s2d = s2d,s1d
                                s1hpop,s2hpop,s1dpop,s2dpop = gff[s1h][self.position],gff[s2h][self.position],gff[s1d][self.position],gff[s2d][self.position]
                                s1hchr,s2hchr,s1dchr,s2dchr = gff[s1h]['chr'],gff[s2h]['chr'],gff[s1d]['chr'],gff[s2d]['chr']         
                                s1hpo,s2hpo,s1dpo,s2dpo = se_dics[num][s1hchr]['start']+s1hpop,se_dics[num+1][s2hchr]['start']+s2hpop,se_dics[num][s1dchr]['start']+s1dpop,se_dics[num+1][s2dchr]['start']+s2dpop
                                s1hpop,s2hpop,s1dpop,s2dpop = s1hpo*step,s2hpo*step,s1dpo*step,s2dpo*step
                                x,y = self.drow_coli(s1hpop,s2hpop,s1dpop,s2dpop,h0,h0-mh,True)
                                plt.fill(x, y, facecolor=chro_color,alpha=self.alpha_f)

                            elif s1hchr in chrlists[num+1] and s2dchr in chrlists[num]:# 共线性第一列就是第二排
                                # s1h,s2h = s2h,s1h
                                # s1d,s2d = s2d,s1d
                                s1hpop,s2hpop,s1dpop,s2dpop = gff[s1h][self.position],gff[s2h][self.position],gff[s1d][self.position],gff[s2d][self.position]
                                s1hchr,s2hchr,s1dchr,s2dchr = gff[s1h]['chr'],gff[s2h]['chr'],gff[s1d]['chr'],gff[s2d]['chr']         
                                s1hpo,s2hpo,s1dpo,s2dpo = se_dics[num+1][s1hchr]['start']+s1hpop,se_dics[num][s2hchr]['start']+s2hpop,se_dics[num+1][s1dchr]['start']+s1dpop,se_dics[num][s2dchr]['start']+s2dpop
                                s1hpop,s2hpop,s1dpop,s2dpop = s1hpo*step,s2hpo*step,s1dpo*step,s2dpo*step
                                x,y = self.drow_coli(s1hpop,s2hpop,s1dpop,s2dpop,h0-mh,h0,False)
                                plt.fill(x, y, facecolor=chro_color,alpha=self.alpha_f)

                        else:
                            if s1hchr in chrlists[num] and s2dchr in chrlists[num+1]:# 共线性第一列就是第二排
                                s1hpop,s2hpop,s1dpop,s2dpop = gff[s1h][self.position],gff[s2h][self.position],gff[s1d][self.position],gff[s2d][self.position]
                                s1hchr,s2hchr,s1dchr,s2dchr = gff[s1h]['chr'],gff[s2h]['chr'],gff[s1d]['chr'],gff[s2d]['chr']         
                                s1hpo,s2hpo,s1dpo,s2dpo = se_dics[num][s1hchr]['start']+s1hpop,se_dics[num+1][s2hchr]['start']+s2hpop,se_dics[num][s1dchr]['start']+s1dpop,se_dics[num+1][s2dchr]['start']+s2dpop
                                s1hpop,s2hpop,s1dpop,s2dpop = s1hpo*step,s2hpo*step,s1dpo*step,s2dpo*step
                                x,y = self.drow_coli(s1hpop,s2hpop,s1dpop,s2dpop,h0,h0-mh,True)
                                plt.fill(x, y, facecolor=chro_color,alpha=self.alpha_f)

                            elif s1hchr in chrlists[num+1] and s2dchr in chrlists[num]:# 共线性第一列就是第二排
                                s1hpop,s2hpop,s1dpop,s2dpop = gff[s1h][self.position],gff[s2h][self.position],gff[s1d][self.position],gff[s2d][self.position]
                                s1hchr,s2hchr,s1dchr,s2dchr = gff[s1h]['chr'],gff[s2h]['chr'],gff[s1d]['chr'],gff[s2d]['chr']         
                                s1hpo,s2hpo,s1dpo,s2dpo = se_dics[num+1][s1hchr]['start']+s1hpop,se_dics[num][s2hchr]['start']+s2hpop,se_dics[num+1][s1dchr]['start']+s1dpop,se_dics[num][s2dchr]['start']+s2dpop
                                s1hpop,s2hpop,s1dpop,s2dpop = s1hpo*step,s2hpo*step,s1dpo*step,s2dpo*step
                                x,y = self.drow_coli(s1hpop,s2hpop,s1dpop,s2dpop,h0-mh,h0,False)
                                plt.fill(x, y, facecolor=chro_color,alpha=self.alpha_f)

                    elif s1hchr == s2hchr and s1hchr in chrlists[num] and s2dchr in chrlists[num+1]:
                        s1hpo,s2hpo,s1dpo,s2dpo = se_dics[num][s1hchr]['start']+s1hpop,se_dics[num][s2hchr]['start']+s2hpop,se_dics[num][s1dchr]['start']+s1dpop,se_dics[num][s2dchr]['start']+s2dpop
                        s1hpop,s2hpop,s1dpop,s2dpop = s1hpo*step,s2hpo*step,s1dpo*step,s2dpo*step
                        x,y = self.drow_coli(s1hpop,s2hpop,s1dpop,s2dpop,h0,h0-mh,False)
                        plt.fill(x, y, facecolor=chro_color,alpha=self.alpha_f)

            else:
                for s1h,s2h in colineartlys[num]:
                    if self.blast_reverse == "True":
                        s1h,s2h = s2h,s1h
                    s1hpop,s2hpop = gff[s1h][self.position],gff[s2h][self.position]
                    s1hchr,s2hchr = gff[s1h]['chr'],gff[s2h]['chr']
                    if s1hchr in chr_color.keys():
                        chro_color = chr_color[s1hchr]
                    elif s2hchr in chr_color.keys():
                        chro_color = chr_color[s2hchr]
                    else:
                        chro_color = self.linecs.split(',')[num]

                    if s1hchr != s2hchr:# 如果当前block在不同染色体上
                        if s1hchr in chrlists[num] and s2hchr in chrlists[num+1]:# 共线性第一列就是第二排
                            s1hpop,s2hpop = gff[s1h][self.position],gff[s2h][self.position]
                            s1hchr,s2hchr = gff[s1h]['chr'],gff[s2h]['chr']
                            s1hpo,s2hpo = se_dics[num][s1hchr]['start']+s1hpop,se_dics[num+1][s2hchr]['start']+s2hpop
                            s1hpop,s2hpop = s1hpo*step,s2hpo*step
                            x,y = self.drow_coli0(s1hpop,s2hpop,h0,h0-mh,True)
                            plt.plot(x, y, color=chro_color,alpha=self.alpha_f,linewidth=self.linewidth)
                        elif s1hchr in chrlists[num+1] and s2hchr in chrlists[num]:# 共线性第一列就是第二排
                            s1hpop,s2hpop = gff[s1h][self.position],gff[s2h][self.position]
                            s1hchr,s2hchr = gff[s1h]['chr'],gff[s2h]['chr']
                            s1hpo,s2hpo = se_dics[num+1][s1hchr]['start']+s1hpop,se_dics[num][s2hchr]['start']+s2hpop
                            s1hpop,s2hpop = s1hpo*step,s2hpo*step
                            x,y = self.drow_coli0(s1hpop,s2hpop,h0-mh,h0,False)
                            plt.plot(x, y, color=chro_color,alpha=self.alpha_f,linewidth=self.linewidth)

            if self.focus_list != 'None':
                for s1,s2 in focus_list2[num]:
                    if focus_pair == []:
                        try:
                            color = focus_listc[s1]
                        except:
                            color = focus_listc[s2]
                        pass
                    else:
                        if (s1,s2) in focus_pair or (s2,s1) in focus_pair:
                            try:
                                color = focus_listc[s1+'lanmeifang'+s2]
                            except:
                                color = focus_listc[s2+'lanmeifang'+s1]
                        else:
                            continue
                    if self.blast_reverse == "True":
                        s1,s2 = s2,s1
                    s1pop,s2pop = gff[s1][self.position],gff[s2][self.position]
                    s1chr,s2chr = gff[s1]['chr'],gff[s2]['chr']
                    if s1chr not in (chrlists[num]+chrlists[num+1]) or s2chr not in (chrlists[num]+chrlists[num+1]):# 排除不在绘制范围内的共线性
                        continue
                    if s1chr != s2chr:# 如果当前block在不同染色体上
                        if s1chr in chrlists[num] and s2chr in chrlists[num+1]:# 共线性第一列就是第二排
                            s1pop,s2pop = gff[s1][self.position],gff[s2][self.position]
                            s1chr,s2chr = gff[s1]['chr'],gff[s2]['chr']       
                            s1po,s2po = se_dics[num][s1chr]['start']+s1pop,se_dics[num+1][s2chr]['start']+s2pop
                            s1pop,s2pop = s1po*step,s2po*step
                            x,y = plot_bezier_curve(s1pop,h0-(self.h)*0.1,s2pop,h0-mh+(self.h)*1.1)
                            plt.plot(x, y, color = color, linewidth = self.fucus_lw, solid_capstyle='round')
                        elif s1chr in chrlists[num+1] and s2chr in chrlists[num]:# 共线性第一列就是第二排
                            s1pop,s2pop = gff[s1][self.position],gff[s2][self.position]
                            s1chr,s2chr = gff[s1]['chr'],gff[s2]['chr']       
                            s1po,s2po = se_dics[num+1][s1chr]['start']+s1pop,se_dics[num][s2chr]['start']+s2pop
                            s1pop,s2pop = s1po*step,s2po*step
                            x,y = plot_bezier_curve(s1pop,h0-mh+(self.h)*1.1,s2pop,h0-(self.h)*0.1)
                            plt.plot(x, y, color = color, linewidth = self.fucus_lw, solid_capstyle='round')
                    elif s1chr == s2chr and s1chr in chrlists[num] and s2chr in chrlists[num+1]:
                        s1po,s2po = se_dics[num][s1chr]['start']+s1pop,se_dics[num][s2chr]['start']+s2pop
                        s1pop,s2pop = s1po*step,s2po*step
                        x,y = plot_bezier_curve(s1pop,h0-(self.h)*0.1,s2pop,h0-mh+(self.h)*1.1)
                        plt.plot(x, y, color = color, linewidth = self.fucus_lw, solid_capstyle='round')
            h0 -= mh


        ax.set_xlim(0, 1)
        # 关闭坐标轴
        ax.axis('off')
        plt.savefig(self.savefile,dpi=1000 , bbox_inches = 'tight')#