# -*- encoding: utf-8 -*-
'''
@File        :Ks_block.py
@Time        :2021/09/28 11:20:59
@Author        :charles kiko
@Version        :1.0
@Contact        :charles_kiko@163.com
@Desc        :block的ks
'''


from scipy.stats.kde import gaussian_kde
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from famCircle.bez import *
import sys
import csv
from matplotlib.font_manager import FontProperties
from matplotlib import rcParams
from famCircle.bez import get_path


class Ks_block():
    def __init__(self, options):
        self.vertical = "False"
        self.model = "NG86"
        self.dpi = 600
        self.peaks = 'None'
        self.bins = 250
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)
        self.bins = int(self.bins)

    def readks(self):
        read_ks = {}
        f = open(self.ks, 'r', encoding='utf-8')
        for row in f:
            if row[0] != '#' and row[0] != '\n':
                row = row.strip('\n').split('\t')
                if row[0] != 'id1' and len(row) != 2:
                    pair = str(row[0]) + '_' + str(row[1])
                    if self.model == 'NG86':
                        lt = float(row[3])
                    elif self.model == 'YN00':
                        lt = float(row[5])
                    else:
                        print('未能解析ks文件')
                        exit()
                    if lt < float(self.area.split(',')[0]) or lt > float(self.area.split(',')[1]):
                        continue
                    read_ks[pair] = lt
        return read_ks

    def readblast(self):
        one_gene = []
        lt00 = []
        alphagenepairs = open(self.genepairs, 'r', encoding='utf-8')
        if self.genepairsfile_type == 'famCircle':
            for row in alphagenepairs:
                if (row[0] == '\n'):
                    continue
                elif (row[0] == '#'):
                    if len(lt00) == 0:
                        pass
                    else:
                        one_gene.append(lt00)
                        lt00 = []
                    lt = row.strip('\n').split(' ')
                    block = {}
                    for i in lt:
                        if '=' in str(i):
                            lt0 = i.split('=')
                            block[str(lt0[0])] = str(lt0[1])
                    N = int(block['N'])
                    if N >= int(self.cutoff):
                        self.class1 = True
                    else :
                        self.class1 = False
                else:
                    if self.class1:
                        lt = row.strip('\n').split()
                        id1, id2 = str(lt[0]),str(lt[1])
                        lt00.append([id1,id2])
        elif self.genepairsfile_type == 'WGDI':
            for row in alphagenepairs:
                if (row[0] == '\n'):
                    continue
                elif (row[0] == '#'):
                    if len(lt00) == 0:
                        pass
                    else:
                        one_gene.append(lt00)
                        # print(lt00)
                        lt00 = []
                    lt = row.strip('\n').split(' ')
                    block = {}
                    for i in lt:
                        if '=' in str(i):
                            lt0 = i.split('=')
                            block[str(lt0[0])] = str(lt0[1])
                    N = int(block['N'])
                    if N >= int(self.cutoff):
                        self.class1 = True
                    else :
                        self.class1 = False
                else:
                    if self.class1:
                        lt = row.strip('\n').split()
                        id1, id2 = str(lt[0]),str(lt[2])
                        lt00.append([id1,id2])
        elif self.genepairsfile_type == 'ColinearScan':
            for row in alphagenepairs:
                if (row[0] == '\n',row[0] == '+'):
                    continue
                elif (row[:3] == 'the'):
                    if len(lt00) == 0:
                        pass
                    else:
                        one_gene.append(lt00)
                        lt00 = []
                    lt = row.strip('\n').split()
                    N = int(lt[-1])
                    if N >= int(self.cutoff):
                        self.class1 = True
                    else :
                        self.class1 = False
                else:
                    if self.class1:
                        lt = row.strip('\n').split()
                        id1, id2 = str(lt[0]),str(lt[2])
                        lt00.append([id1,id2])
        elif self.genepairsfile_type == 'MCScanX':
            for row in alphagenepairs:
                if (row[0] == '\n'):
                    continue
                elif(row[:12] == '## Alignment'):
                    if len(lt00) == 0:
                        pass
                    else:
                        one_gene.append(lt00)
                        lt00 = []
                    lt = row.strip('\n').split(' ')
                    block = {}
                    for i in lt:
                        if '=' in str(i):
                            lt0 = i.split('=')
                            block[str(lt0[0])] = str(lt0[1])
                    N = int(block['N'])
                    if N >= int(self.cutoff):
                        self.class1 = True
                    else :
                        self.class1 = False
                elif ('#' not in row):
                    if self.class1:
                        lt = row.strip('\n').split()
                        # print(lt)
                        if len(lt) == 5:
                            id1, id2 = str(lt[2]),str(lt[3])
                        elif len(lt) == 4:
                            id1, id2 = str(lt[1]),str(lt[2])
                        elif len(lt) == 6:
                            id1, id2 = str(lt[3]),str(lt[4])
                        else:
                            # print(row)
                            print('Parse error!')
                            exit()
                        lt00.append([id1,id2])
        elif self.genepairsfile_type == 'BLAST':
            alphagenepairs = pd.read_csv(self.genepairs,header = None, sep='\t')
            alphagenepairs = alphagenepairs.drop(alphagenepairs[alphagenepairs[0] == alphagenepairs[1]].index)
            dic0 = alphagenepairs.groupby(0).groups
            job = len(dic0)
            i = 0
            for key in dic0.keys():
                p = int(((i+1)/job)*30)
                print("\r["+"*"*p+" "*(30-p)+"] "+str(i)+"/"+str(job)+" Data cleansing ",end="")
                local = alphagenepairs.loc[dic0[key]].sort_values(by=[0,3],ascending= [True, False])
                local.reset_index(drop=True,inplace=True)
                local = local.head(int(self.cutoff))
                # local[13] = [local[0],local[1]]
                local[13] = local.apply(lambda row: [row[0], row[1]], axis=1)
                lt = local[13].to_list()
                one_gene.append(lt)
                i += 1
        else:
            print('genepairsfile_type error: File Format not recognized!')
            exit()
        if self.genepairsfile_type == 'BLAST' or len(lt0) == 0:
            pass
        else:
            one_gene.append(lt00)
            lt00 = []
        return one_gene

        # 绘制核密度曲线图


    def KdePlot(self,list0,font_prop):
        vertical = 'False'
        # y = list0[0]# 平均数
        x = list0# 中位数x = list(map(int, x))
        fig = plt.figure(figsize=(20,10),dpi=800)


        plt.grid(c='grey',ls='--',linewidth=0.3)
        plt.hist(x, bins=self.bins, density=True, alpha=0.4,color = self.fillc)
        dist_space = np.linspace(float(self.area.split(',')[0]),float(self.area.split(',')[1]), self.bins)
        kdemedian = gaussian_kde(x)
        kdemedian.set_bandwidth(bw_method=kdemedian.factor / 3.)
        xx = list(dist_space)[list(kdemedian(dist_space)).index(max(list(kdemedian(dist_space))))]
        xx0 = round(xx,4)
        plt.plot(dist_space, kdemedian(dist_space), color=self.linec,
                 label=self.ks_value+'_KDE_' + str(xx0))

        plt.legend(prop=font_prop,fontsize=24)
        plt.title('Kernel Density Plot of Ks Values', fontproperties=font_prop,fontsize=28)# 设置图片标题
        plt.xlabel('ks', fontproperties=font_prop,fontsize=26)# 设置 x 轴标签
        plt.ylabel('density', fontproperties=font_prop,fontsize=26)# 设置 y 轴标签

        # 设置坐标轴刻度的字体大小
        plt.tick_params(axis='both', which='major', labelsize=20)  # 设置主刻度字体大小
        plt.tick_params(axis='both', which='minor', labelsize=20)  # 设置次刻度字体大小
        # # 为刻度标签应用字体属性
        # for tick in plt.get_xticklabels() + plt.get_yticklabels():
        #     tick.set_fontproperties(font_prop)

    # def multi_norm(self,x, *params):
    #     result = 0
    #     for i in range(0, len(params), 3): 
    #         mean, std_dev, weight = params[i], params[i+1], params[i+2] 
    #         result += mean * np.exp(-(x - std_dev)**2 / (2 * weight**2))
    #     return result

    def multi_norm(self, x, *params):
        n = len(params) // 3
        result = np.zeros_like(x)
        for i in range(n):
            mean = params[3*i]
            std = params[3*i + 1]
            amplitude = params[3*i + 2]
            result += amplitude * np.exp(-((x - mean) ** 2) / (2 * std ** 2))
        return result


    def AMPD(self,data):
        """
        实现AMPD算法
        :param data: 1-D numpy.ndarray 
        :return: 波峰所在索引值的列表
        """ 
        data = np.array(data)
        p_data = np.zeros_like(data, dtype=np.int32)
        count = data.shape[0]
        arr_rowsum = []
        for k in range(1, count // 2 + 1):
            row_sum = 0
            for i in range(k, count - k):
                if data[i] > data[i - k] and data[i] > data[i + k]:
                    row_sum -= 1
            arr_rowsum.append(row_sum)
        min_index = np.argmin(arr_rowsum)
        max_window_length = min_index
        for k in range(1, max_window_length + 1):
            for i in range(k, count - k):
                if data[i] > data[i - k] and data[i] > data[i + k]:
                    p_data[i] += 1
        return np.where(p_data == max_window_length)[0]

    def Histplot(self,ks_data,font_prop):
        plt.figure(figsize=(20,10),dpi=800)
        # 获取字体名称
        font_name = font_prop.get_name()
        # 设置全局字体属性
        rcParams['font.family'] = font_name
        rcParams['font.size'] = 12
        rcParams['xtick.labelsize'] = 20
        rcParams['ytick.labelsize'] = 20

        plt.grid(c='grey',ls='--',linewidth=0.3)
        from scipy.optimize import curve_fit
        from scipy.signal import argrelextrema
        plt.hist(ks_data, bins=self.bins, density=True, alpha=0.4,color = self.fillc)
        band = float(self.area.split(',')[1])/20

        if self.peaks == 'None':
            hist, bin_edges = np.histogram(ks_data, bins=int(self.bins/2), density=True)
            px = self.AMPD(hist)
            peaks = [px]
            data = ks_data
            # guess = []
            # for i in peaks[0]:
            #     guess = guess + [bin_edges[i],0.3,hist[i]]

            # hist, bin_edges = np.histogram(ks_data, bins=int(self.bins/2), density=True)
            # from scipy.signal import find_peaks

            # # 示例数组，包含多个峰值
            # arr = np.array(hist)

            # # 使用 find_peaks 查找峰值的位置
            # peaks, _ = find_peaks(arr)
            # if len(peaks) > 2:
            #     hist0 = [hist[i] for i in peaks]
            #     # 示例数组，包含多个峰值
            #     arr = np.array(hist0)
            #     # 使用 find_peaks 查找峰值的位置
            #     peaks, _ = find_peaks(arr)
            #     if len(peaks) > 3:
            #         peaks = [np.where(hist == max(hist))[0]]
            # lt = [(bin_edges[i]+bin_edges[i+1])/3 for i in peaks]
            # data = []
            # for i in lt:
            #     data += [x for x in ks_data if i-band <= x <= i+band]
            # # print(len(data),len(ks_data))
            # data0 = []
            # for i in range(round(len(ks_data)/len(data))):
            #     data0 += data 
            # data = data0
            # for i in range(int((float(self.area.split(',')[1])-float(self.area.split(',')[0]))/0.001)):
            #     data.append(i*0.001)
            # print(len(data))
        else:
            lt = [float(i) for i in self.peaks.split(',')]
            data = []

            for i in lt:
                data += [x for x in ks_data if i-band <= x <= i+band]
            # print(len(data),len(ks_data))
            data0 = []
            for i in range(round(len(ks_data)/len(data))):
                data0 += data 
            data = data0
            for i in range(int((float(self.area.split(',')[1])-float(self.area.split(',')[0]))/0.001)):
                data.append(i*0.001)
            # print(len(data))

        hist, bin_edges = np.histogram(data, bins=self.bins, density=True)
        # print(hist,bin_edges)
        px = self.AMPD(hist)
        # print(px)
        peaks = [px]
        guess = []
        for i in peaks[0]:
            guess = guess + [bin_edges[i],0.3,hist[i]]
        # print(guess)
        # exit()
        max_peaks = 4  # 设置最多使用的峰值数量
        guess = guess[:3 * max_peaks]

        # 创建直方图的bins的中心点
        bin_edges_center = (bin_edges[:-1] + bin_edges[1:]) / 2

        n = len(hist)
        maxfev = int(np.ceil(100000 * n))
        params, cov = curve_fit(self.multi_norm, bin_edges_center, hist, guess, maxfev=maxfev)

        plt.plot(bin_edges_center, self.multi_norm(bin_edges_center, *params), label='AMPD fit Gaussian distribution',color = self.linec) #, "b"+marker[i]+":"
        

        plt.legend(prop=font_prop,fontsize=24)
        plt.title('Histogram of Ks Values with Normal Distribution Fit', fontproperties=font_prop,fontsize=28)# 设置图片标题
        plt.xlabel('ks', fontproperties=font_prop,fontsize=26)# 设置 x 轴标签
        plt.ylabel('density', fontproperties=font_prop,fontsize=26)# 设置 y 轴标签

        # # 设置坐标轴刻度的字体大小
        # plt.tick_params(axis='both', which='major', labelsize=20)  # 设置主刻度字体大小
        # plt.tick_params(axis='both', which='minor', labelsize=20)  # 设置次刻度字体大小
        # # 为刻度标签应用字体属性
        # for tick in plt.get_xticklabels() + plt.get_yticklabels():
        #     tick.set_fontproperties(font_prop)

    def run(self):
        read_ks = self.readks()
        blast = self.readblast()
        path = get_path()
        font_path = os.path.join(path, 'example/arial.ttf')
        font_prop = FontProperties(fname=font_path)
        # 获取字体名称
        font_name = font_prop.get_name()
        # 设置全局字体属性
        rcParams['font.family'] = font_name
        rcParams['font.size'] = 12
        rcParams['xtick.labelsize'] = 20
        rcParams['ytick.labelsize'] = 20

        # # 确认全局设置
        # print("rcParams 设置：")
        # print(rcParams['font.family'])
        # print(rcParams['font.size'])
        # print(rcParams['xtick.labelsize'])
        # print(rcParams['ytick.labelsize'])
        average_list = []
        median_list = []
        for i in blast:
            lt = []
            for j in i:
                if str(j[0]) + '_' + str(j[1]) in read_ks.keys():
                    ks_v = read_ks[str(j[0]) + '_' + str(j[1])]
                elif str(j[1]) + '_' + str(j[0]) in read_ks.keys():
                    ks_v = read_ks[str(j[1]) + '_' + str(j[0])]
                else:
                    continue
                lt.append(ks_v)
            if len(lt) == 0:
                continue
            average_list.append(np.mean(lt))
            median_list.append(np.median(lt))
        if self.ks_value == 'median':
            ks_list = median_list
        elif self.ks_value == 'mean':
            ks_list = average_list
        else:
            print('Please set the parameter ks_value = median/mean, typically chosen as median.')
        if self.plot_type == 'Kernel':
            # print('我在这里呀')
            self.KdePlot(ks_list,font_prop)
        elif self.plot_type == 'Gaussian':
            self.Histplot(ks_list,font_prop)
        plt.savefig(self.savefile, dpi = int(self.dpi), bbox_inches = 'tight')# 存储图片
        sys.exit(0)