### draw chromosome homology famCircle for eudicots with VV as reference, if VV is provided. If no vv chro is involved, it draws famCircle for other species
# -*- coding: UTF-8 -*-


import re
import sys
from math import *

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import *
from matplotlib.patches import Circle, Ellipse
from pylab import *
from scipy.optimize import curve_fit
import math
import matplotlib.mlab as mlab
from scipy.stats import norm

from famCircle.bez import *
matplotlib.rcParams['font.sans-serif'] = ['KaiTi']

class lookKs():
    def __init__(self, options):
    ##### circle parameters
        self.y1 = []
        self.y2 = []
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)


    def ks_scatter_plot(self, x, y, z):
        z = np.array(z) # 这里填入你的数据list 如果已经是array格式就不用转化了
        z = list(map(float, z))
        mu =np.mean(z) #计算均值 
        sigma =np.std(z) 
        num_bins = len(x) #直方图柱子的数量 
        n, bins, patches = plt.hist(z, num_bins, log = False, density=1, alpha=0.75) 
        #直方图函数，x为x轴的值，normed=1表示为概率密度，即和为一，绿色方块，色深参数0.5.返回n个概率，直方块左边线的x值，及各个方块对象 
        y1 = norm.pdf(bins, mu, sigma)#拟合一条最佳正态分布曲线y 
        plt.grid(True)
        plt.plot(bins, y1, 'r--') #绘制y的曲线 

    def run(self):
        fpgenefamilyinf = open(self.ks, 'r', encoding='utf-8')
        for row in fpgenefamilyinf:
            if (row.split()[0] == 'id1'):
                continue
            gene1, gene2, Ka, Ks = row.split()[0],row.split()[1],row.split()[2],row.split()[3]
            ka = '%.2f' % eval(Ka)
            ks = '%.2f' % eval(Ks)
            self.y1.append(ks)
            self.y2.append(ka)
        fpgenefamilyinf.close()
        MAX = float(max(self.y1))
        MIN = float(min(self.y1))
        set1 = set(self.y1)
        dict1 = {}
        for item in set1:
            dict1.update({item:self.y1.count(item)})
        set2 = set(self.y2)
        dict2 = {}
        for item in set2:
            dict2.update({item:self.y2.count(item)})
        # print(dict1,dict2)
        x = sorted(dict1.keys())
        y = []
        for i in x:
            y.append(dict1[i])
        # print(x, y)
        plt.figure() #初始化一张图 
        plt.xlabel('Ks')  
        plt.ylabel('frequency')  
        plt.title(r'Ks distribution map') #+citys[i])  
        plt.xticks(np.arange(MIN,MAX,0.5),rotation=25)
        self.ks_scatter_plot(x, y, self.y1)
        savefig(self.savefile, dpi=500)
        sys.exit(0)
