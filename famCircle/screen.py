# -*- encoding: utf-8 -*-
'''
@File        :screen.py
@Time        :2021/09/28 09:03:54
@Author        :charles kiko
@Version        :1.0
@Contact        :charles_kiko@163.com
@Desc        :分染色体展示家族分布折线图
'''


from famCircle.bez import *
import os
import sys
from pylab import *

class screen():
    def __init__(self, options):
        self.marker = ['.','o','^','v','<','>','s','+','x']
        self.color = ['red', 'fuchsia', 'aquamarine', 'orangered', 'violet', 'lime', 'chocolate', 'blueviolet', 'limegreen', 'orange', 'royalblue', 'gold', 'dodgerblue', 'yellow', 'cyan', 'yellowgreen', 'teal', 'palegreen']
        self.ls = ['-','-','-']
        self.series = "25"
        self.position = 'order'
        self.dpi = 600
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)

    def run_plot(self,x_axis_data,y_axis_data,name,tup,domain):
        linestyle = str(tup[0])
        marker = str(tup[1])
        color = str(tup[2])
        # plot中参数的含义分别是横轴值，纵轴值，线的形状，颜色，透明度,线的宽度和标签
        plt.plot(x_axis_data, y_axis_data, marker = marker, markersize = 0.5, linestyle = linestyle, color = color, alpha = 0.5, linewidth = 1, label = domain)
        # 显示标签，如果不加这句，即使在plot中加了label='一些数字'的参数，最终还是不会显示标签
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
        plt.xlabel('chr')
        plt.ylabel('gene')

    def readgf(self,file):
        f = open(self.domainpath + '/' + file, 'r', encoding='utf-8')
        genef = []
        for row in f:
            row = row.strip()
            if "#" not in row and row != '\n':
                row = row.strip('\n').split()
                id0 = row[0]
                if id0 not in genef:
                    genef.append(id0)
        for i in genef:
            if i == None:
                genef.remove(None)
        return genef

    def numb(self,i):
        ls = self.ls[i % 3]
        marker = self.marker[i % 9]
        color = self.color[i % 18]
        return ls,marker,color

    def run(self):
        if os.path.isdir(self.outpath):
            pass
        else:
            d = os.mkdir(self.outpath)
        gff = read_gff(self.gff)
        lens,chr_list = read_lens(self.lens)
        name = os.listdir(self.domainpath)
        average = sum([lens[i]['end'] for i in lens.keys()])/sum([lens[i]['order'] for i in lens.keys()])
        domains = {}
        for i in range(len(name)):
            domain = name[i]
            nm = str(str(domain).split('.')[0])
            domains[nm] = self.readgf(domain) + [self.numb(i)]
        data = {}
        for key in domains.keys():# 遍历结构域
            chrx = {}
            sx = domains[key][-1]
            lt = domains[key][:-1]
            chrx['att'] = sx
            for j in range(len(lt) - 1):# 遍历基因
                jishu = {}
                i = lt[j]
                chrg = gff[i]['chr']
                start = gff[i]['end']
                jishug = int(int(start)/(int(self.series) * average))
                # print(chrg,jishug)
                if chrg not in chrx.keys():
                    jishu[jishug] = 1
                    chrx[chrg] = jishu
                else:
                    if jishug in chrx[chrg].keys():
                        chrx[chrg][jishug] += 1
                    else:
                        chrx[chrg][jishug] = 1
            for i in list(chrx.keys()):
                try:
                    length0 = max(list(chrx[i].keys()))
                    # print(length0)
                    for j in range(length0 + 1):
                        if j in list(chrx[i].keys()):
                            pass
                        else:
                            chrx[i][j] = 0
                except:
                    pass
            data[key] = chrx
            chrx = {}
        chry = []
        for domain in data.keys():
            for ch in data[domain].keys():
                if str(ch) == 'att':
                    continue
                else:
                    chry.append(ch)
            break
        domains = list(data.keys())
        for i in chry:
            if i not in chr_list:
                continue
            for j in domains:
                tup = data[j]['att']
                name = str(i) + '_' + str(j) + '.png'
                x = sorted(data[j][i].keys())
                # x = list(data[domain][ch].keys())
                y = []
                for o in x:
                    y.append(data[j][i][o])
                self.run_plot(x,y,name,tup,j)
            plt.title(i)
            plt.savefig(self.outpath + '/' +str(i), dpi = int(self.dpi),bbox_inches='tight')  # 保存该图片
            plt.cla()
        sys.exit(0)
