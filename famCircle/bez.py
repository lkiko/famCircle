# -*- encoding: utf-8 -*-
'''
@File        :bez.py
@Time        :2021/09/28 11:26:28
@Author        :charles kiko
@Version        :1.0
@Contact        :charles_kiko@163.com
@Desc        :None
'''

import configparser
import subprocess
import os
import re
import math
from math import *
import numpy as np
import pandas as pd
import famCircle
from io import StringIO
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
# from Bio.Align.Applications import MafftCommandline, MuscleCommandline, ClustalwCommandline
import codecs
from tqdm import trange
import gc
import matplotlib.pyplot as plt
from matplotlib.patches import *
from matplotlib.patches import Circle, Ellipse
from pylab import *
from collections import Counter
from matplotlib.font_manager import FontProperties

# Bezier functions
#2022.8.15 软件维护

def read_lens(lens):
    chr_dic = {}
    chr_list = []
    lens = pd.read_csv(lens,header = None, sep='\t', comment='#')
    lens[0] = lens[0].astype(str)
    chr_list = lens[0].to_list()
    for index,row in lens.iterrows():
        dic0 = {}
        dic0['end'] = row[1]
        try:
            dic0['order'] = row[2]
        except:
            pass
        # dic0['order'] = row[2]
        try:
            dic0['rev'] = row[3]
        except:
            pass
        chr_dic[row[0]] = dic0
    return chr_dic,chr_list

def get_path():
    from pathlib import Path
    # 获取当前脚本的路径
    script_path = Path(__file__).resolve()
    # 获取当前脚本所在的目录
    script_dir = script_path.parent
    # print(f"Script path: {script_path}")
    # print(f"Script directory: {script_dir}")
    return script_dir

def read_gff(gff,chr_list):
    gene = {}
    gff = pd.read_csv(gff,header = None, sep='\t', comment='#')
    gff[0] = gff[0].astype(str)
    # print(gff)
    for index,row in gff.iterrows():
        if row[0] not in chr_list:
            continue
        dic0 = {}
        dic0['chr'] = (row[0])
        dic0['end'] = row[2]
        dic0['order'] = row[5]
        gene[row[1]] = dic0
    return gene

def read_family1(file):
    family = []
    family_dic = {}
    dic_motif = {}
    fam = pd.read_csv(file,header = None, sep='\t')
    fam.drop_duplicates(subset=[0], keep='first', inplace=True)# 去重
    family = fam[0].to_list()
    for index,row in fam.iterrows():
        if fam.shape[1] == 1:
            family_dic[row[0]] = '#4d5aaf'
            dic_motif['None'] = '#4d5aaf'
        else:
            family_dic[row[0]] = row[1]
            if row[2] not in dic_motif.keys():
                dic_motif[row[2]] = row[1]
    return family,family_dic,dic_motif

def drop_none(fam,gff):
    fam_l = []
    for i in fam:# 去除不在染色体上的家族成员
        if i not in gff.keys():
            continue
        else:
            fam_l.append(i)
    return fam_l

def read_ks1(file0,file_type,family):
    if file_type == 'ks':
        ks_f = pd.read_csv(file0,header = 0, sep='\t', comment='#')
        ks_f = ks_f.drop(ks_f[(ks_f["id1"].apply(lambda x: x not in family)) | (ks_f["id2"].apply(lambda x: x not in family))].index)# 条件过滤
        # print(a)
        # print(ks_f)
        ks_f["all"] = ks_f["id1"].map(str)+"%"+ks_f["id2"].map(str)
        ks_f["all_"] = ks_f["all"].apply(scale)
        ks_df = ks_f.drop_duplicates("all_",keep='first')
        ks_ = ks_df[["all","ks_NG86"]]# 保留两列
        res = ks_.sort_values(by="ks_NG86", ascending=False)#降序
        # print(res)
        ks_l = res["all"].to_list()
        ks = res.set_index("all")["ks_NG86"].to_dict()
    elif file_type == 'blast':
        ks_f = pd.read_csv(file0,header = None, sep='\t')
        ks_f = ks_f.drop(ks_f[(ks_f[0].apply(lambda x: x not in family)) | (ks_f[1].apply(lambda x: x not in family))].index)# 条件过滤
        # print(a)
        ks_f[12] = ks_f[0].map(str)+"%"+ks_f[1].map(str)
        ks_f[13] = ks_f[12].apply(scale)
        ks_df = ks_f.drop_duplicates(13,keep='first')
        ks_ = ks_df[[12,2]]# 保留两列
        res = ks_.sort_values(by=2, ascending=True)#升序
        # print(res)
        ks_l = res[12].to_list()
        ks = res.set_index(12)[2].to_dict()
    return ks,ks_l

def read_coliearity_sline1(file,N,genepairsfile_type,chr_list):
    # block = {}
    data = []
    alphagenepairs = open(file, 'r', encoding='utf-8').read()
    if genepairsfile_type == 'WGDI':
        file = alphagenepairs.split('# Alignment ')[1:]
        n = 0
        for local in file:
            local = local.strip('\n').split('\n')
            head = local[0].split()
            chro = head[-2].split('&')
            if set(chro).issubset(set(chr_list)):
                n += 1
                local = local[1:]
                data_lt = []
                f_data = [local[0].split()[0],local[0].split()[2],local[-1].split()[0],local[-1].split()[2]]
                if len(local) < N:
                    continue
                for line in local:
                    line = line.split()
                    data.append([line[0],line[2]])
    return data


def read_coliearity_sline(file,N,genepairsfile_type,chr_list):
    block = {}
    data = {}
    alphagenepairs = open(file, 'r', encoding='utf-8').read()
    if genepairsfile_type == 'WGDI':
        file = alphagenepairs.split('# Alignment ')[1:]
        n = 0
        for local in file:
            local = local.strip('\n').split('\n')
            head = local[0].split()
            chro = head[-2].split('&')
            if set(chro).issubset(set(chr_list)):
                n += 1
                local = local[1:]
                data_lt = []
                f_data = [local[0].split()[0],local[0].split()[2],local[-1].split()[0],local[-1].split()[2]]
                if len(local) < N:
                    continue
                for line in local:
                    line = line.split()
                    data_lt.append([line[0],line[2]])
                block[n] = f_data
                data[n] = data_lt
    return data,block


import time
from tqdm import tqdm
def readblast_large(genepairs_,blast_reverse,score_levels_,position_,levels_,tandem_,gff,lens,chr1_start,chr2_start,step1,step2):
    print("Count the number of lines in the BLAST file")
    start_time = time.time()  # 记录开始时间
    total_lines = sum(1 for _ in open(genepairs_,'r'))
    end_time = time.time()  # 记录结束时间
    print(f"Function execution time: {end_time - start_time:.2f} seconds")
    x1,x2,x3,y1,y2,y3 = [],[],[],[],[],[]
    # 初始化进度条
    
    with tqdm(total=total_lines, desc="Reading BLAST", unit="lines") as pbar:
        # 分块读取文件
        chunk_iter = pd.read_csv(genepairs_,header = None, sep='\t', comment='#', chunksize=100000)
        processed_lines = 0
        tail = pd.DataFrame()
        for chunk in chunk_iter:
            chunk.loc[:, 12] = chunk.apply(lambda x : cat_gene(x[0],x[1]), axis=1)
            chunk = chunk.drop_duplicates(subset=12 , keep='first')
            del chunk[12]
            if blast_reverse == 'True':
                chunk[[0, 1]] = chunk[[1, 0]]
            chunk = chunk.drop(chunk[(chunk[0] == chunk[1])].index)#.sort_values(by=[0,11],ascending= [True,False])
            last_key = list(chunk[0].to_list())[-1]
            if not tail.empty:
                chunk = pd.concat([tail, chunk], ignore_index=True)
                chunk.reset_index(drop=True,inplace=True)
            dic = chunk.groupby(0).groups# 按照第几列分组
            
            # 判断是否为最后一块
            if processed_lines == total_lines:
                print("This is the last chunk.")
                tail = pd.DataFrame()
            else:
                tail = chunk.loc[dic[last_key]]
                tail.reset_index(drop=True,inplace=True)
                del dic[last_key]

            for key in dic.keys():
                # print(" ******** ",key," ******** ")
                local = chunk.loc[dic[key]]
                if len(local) > 20:
                    continue
                local.reset_index(drop=True,inplace=True)

                if key not in gff.keys():
                    continue
                if score_levels_ == 'False':
                    if tandem_ == 'True':
                        object0 = local[1].to_list()
                        object0 = list(dict.fromkeys(object0))# 去除tandem保留顺序
                        for i in range(levels_[0]):
                            if i < len(object0):
                                if object0[i] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[object0[i]]['chr'] not in lens.keys():
                                    continue
                                if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                    x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                                else:
                                    x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                                if 'rev' in lens[gff[object0[i]]['chr']].keys() and lens[gff[object0[i]]['chr']]['rev'] == '-':
                                    y = -(chr2_start[gff[object0[i]]['chr']] + lens[gff[object0[i]]['chr']][position_] - gff[object0[i]][position_] + 1)*step2
                                else:
                                    y = -(chr2_start[gff[object0[i]]['chr']] + gff[object0[i]][position_])*step2
                                x1.append(x)
                                y1.append(y)
                        for i in range(levels_[0],levels_[0]+levels_[1]):
                            if i < len(object0):
                                if object0[i] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[object0[i]]['chr'] not in lens.keys():
                                    continue
                                if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                    x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                                else:
                                    x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                                if 'rev' in lens[gff[object0[i]]['chr']].keys() and lens[gff[object0[i]]['chr']]['rev'] == '-':
                                    y = -(chr2_start[gff[object0[i]]['chr']] + lens[gff[object0[i]]['chr']][position_] - gff[object0[i]][position_] + 1)*step2
                                else:
                                    y = -(chr2_start[gff[object0[i]]['chr']] + gff[object0[i]][position_])*step2
                                x2.append(x)
                                y2.append(y)
                        for i in range(levels_[0]+levels_[1],levels_[0]+levels_[1]+levels_[2]):
                            if i < len(object0):
                                if object0[i] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[object0[i]]['chr'] not in lens.keys():
                                    continue
                                if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                    x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                                else:
                                    x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                                if 'rev' in lens[gff[object0[i]]['chr']].keys() and lens[gff[object0[i]]['chr']]['rev'] == '-':
                                    y = -(chr2_start[gff[object0[i]]['chr']] + lens[gff[object0[i]]['chr']][position_] - gff[object0[i]][position_] + 1)*step2
                                else:
                                    y = -(chr2_start[gff[object0[i]]['chr']] + gff[object0[i]][position_])*step2
                                x3.append(x)
                                y3.append(y)
                    else:
                        object0 = local[1].to_list()
                        for i in range(levels_[0]):
                            if i < len(object0):
                                if object0[i] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[object0[i]]['chr'] not in lens.keys():
                                    continue
                                if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                    x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                                else:
                                    x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                                if 'rev' in lens[gff[object0[i]]['chr']].keys() and lens[gff[object0[i]]['chr']]['rev'] == '-':
                                    y = -(chr2_start[gff[object0[i]]['chr']] + lens[gff[object0[i]]['chr']][position_] - gff[object0[i]][position_] + 1)*step2
                                else:
                                    y = -(chr2_start[gff[object0[i]]['chr']] + gff[object0[i]][position_])*step2
                                x1.append(x)
                                y1.append(y)
                        for i in range(levels_[0],levels_[0]+levels_[1]):
                            if i < len(object0):
                                if object0[i] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[object0[i]]['chr'] not in lens.keys():
                                    continue
                                if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                    x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                                else:
                                    x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                                if 'rev' in lens[gff[object0[i]]['chr']].keys() and lens[gff[object0[i]]['chr']]['rev'] == '-':
                                    y = -(chr2_start[gff[object0[i]]['chr']] + lens[gff[object0[i]]['chr']][position_] - gff[object0[i]][position_] + 1)*step2
                                else:
                                    y = -(chr2_start[gff[object0[i]]['chr']] + gff[object0[i]][position_])*step2
                                x2.append(x)
                                y2.append(y)
                        for i in range(levels_[0]+levels_[1],levels_[0]+levels_[1]+levels_[2]):
                            if i < len(object0):
                                if object0[i] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[object0[i]]['chr'] not in lens.keys():
                                    continue
                                if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                    x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                                else:
                                    x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                                if 'rev' in lens[gff[object0[i]]['chr']].keys() and lens[gff[object0[i]]['chr']]['rev'] == '-':
                                    y = -(chr2_start[gff[object0[i]]['chr']] + lens[gff[object0[i]]['chr']][position_] - gff[object0[i]][position_] + 1)*step2
                                else:
                                    y = -(chr2_start[gff[object0[i]]['chr']] + gff[object0[i]][position_])*step2
                                x3.append(x)
                                y3.append(y)
                else:
                    score_levels = [float(i) for i in score_levels_.split(',')]
                    for index,row in local.iterrows():
                        if row[2] >= score_levels[0]:
                            if row[1] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[row[1]]['chr'] not in lens.keys():
                                continue
                            if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                            else:
                                x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                            if 'rev' in lens[gff[row[1]]['chr']].keys() and lens[gff[row[1]]['chr']]['rev'] == '-':
                                y = -(chr2_start[gff[row[1]]['chr']] + lens[gff[row[1]]['chr']][position_] - gff[row[1]][position_] + 1)*step2
                            else:
                                y = -(chr2_start[gff[row[1]]['chr']] + gff[row[1]][position_])*step2
                            x1.append(x)
                            y1.append(y)
                        elif row[2] >= score_levels[1]:
                            if row[1] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[row[1]]['chr'] not in lens.keys():
                                continue
                            if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                            else:
                                x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                            if 'rev' in lens[gff[row[1]]['chr']].keys() and lens[gff[row[1]]['chr']]['rev'] == '-':
                                y = -(chr2_start[gff[row[1]]['chr']] + lens[gff[row[1]]['chr']][position_] - gff[row[1]][position_] + 1)*step2
                            else:
                                y = -(chr2_start[gff[row[1]]['chr']] + gff[row[1]][position_])*step2
                            x2.append(x)
                            y2.append(y)
                        elif row[2] >= score_levels[2]:
                            if row[1] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[row[1]]['chr'] not in lens.keys():
                                continue
                            if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                            else:
                                x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                            if 'rev' in lens[gff[row[1]]['chr']].keys() and lens[gff[row[1]]['chr']]['rev'] == '-':
                                y = -(chr2_start[gff[row[1]]['chr']] + lens[gff[row[1]]['chr']][position_] - gff[row[1]][position_] + 1)*step2
                            else:
                                y = -(chr2_start[gff[row[1]]['chr']] + gff[row[1]][position_])*step2
                            x3.append(x)
                            y3.append(y)
            pbar.update(len(chunk))
            processed_lines += len(chunk)
    return x1,x2,x3,y1,y2,y3


def clusterx_pro(index,j,lt,gene_average,gff,dic0,cluster,magnification):
    for i in range(1,int(cluster)+1):
        if j+i in lt or j-i in lt:
            dow = j-i
            up = j+i
            if dow in lt and abs(gff[dic0[j]]['end']-gff[dic0[j-i]]['end']) < int(cluster)*gene_average*int(magnification):
                return True
            elif up in lt and abs(gff[dic0[j]]['end']-gff[dic0[j+i]]['end']) < int(cluster)*gene_average*int(magnification):# 平均基因长度
                if index != 0 and abs(gff[dic0[j]]['end']-gff[dic0[lt[index-1]]]['end']) < int(cluster)*gene_average*int(magnification):
                    return True
                else: 
                    return False
            else:
                False
    return False

def plot_block_s(start, stop,hight, radius):
    # start, stop measured in radian 染色体
    #print start, stop
    t = np.arange(start, stop, math.pi/720.)
    x_s, y_s = radius*np.cos(t), radius*np.sin(t)
    x_b, y_b = (radius+hight)*np.cos(t), (radius+hight)*np.sin(t)
    # print('***********************************************************')
    # print(x_s,len(x_s))
    # print(x_b,len(x_b))
    x_ = np.append(x_s,x_b[::-1])
    y_ = np.append(y_s,y_b[::-1])
    # print(x_,len(x_))
    x = np.append(x_,x_s[0])
    y = np.append(y_,y_s[0])
    # print(x,len(x))
    # print('***********************************************************')
    return x,y

def plot_block_s_line(start, stop,hight, y):
    # start, stop measured in radian 染色体
    #print start, stop
    x = [start,stop,stop,start,start]
    y = [y,y,y+hight,y+hight,y]
    # print('***********************************************************')
    # print(x,len(x))
    # print('***********************************************************')
    return x,y

#############################绘制染色体步骤


def to_radian(bp, total):
    # from basepair return as radian
    # print ("to_radian", bp, total)
    return np.radians(bp*360./total)

def plot_chr(start,stop,radius,size):
    t = np.arange(start, stop, math.pi/720.)
    x, y = radius*np.cos(t), radius*np.sin(t)

def plot_rec(start, stop, radius,size):
    # start, stop measured in radian 染色体
    t = np.arange(start, stop, math.pi/720.)
    x, y = radius*np.cos(t), radius*np.sin(t)
    x,y = list(x[::-1]),list(y[::-1])

    t = np.arange(start, stop, math.pi/720.)
    x2, y2 = (radius+size)*np.cos(t), (radius+size)*np.sin(t)
    x2,y2 = list(x2),list(y2)
    x = x + x2
    y = y + y2

    x.append(x[0])
    y.append(y[0])

    return x,y

def plot_arc(start, stop, radius,size):
    # start, stop measured in radian 染色体
    #print start, stop
    # print(start,stop)
    # exit()
    t = np.arange(start, stop, math.pi/720.)
    x, y = radius*np.cos(t), radius*np.sin(t)
    x,y = list(x[::-1]),list(y[::-1])
   
    t = np.arange(start, start-math.pi, -math.pi/30.)
    x1, y1 = (size/2)*np.cos(t), (size/2)*np.sin(t)
    middle_r = radius+(size)/2
    x1, y1 = x1 + middle_r*np.cos(start), y1 + middle_r*np.sin(start)
    x1,y1 = list(x1[::-1]),list(y1[::-1])
    x = x + x1
    y = y + y1

    t = np.arange(start, stop, math.pi/720.)
    x2, y2 = (radius+size)*np.cos(t), (radius+size)*np.sin(t)
    x2,y2 = list(x2),list(y2)
    x = x + x2
    y = y + y2

    t = np.arange(stop, stop+math.pi, math.pi/30.)
    x3, y3 = (size/2)*np.cos(t), (size/2)*np.sin(t)
    middle_r = radius+(size)/2
    x3, y3 = x3 + middle_r*np.cos(stop), y3 + middle_r*np.sin(stop)
    x3,y3 = list(x3),list(y3)
    x,y = x + x3,y + y3
    x.append(x[0])
    y.append(y[0])

    return x,y

def line_chr1(h,xmin,xmax,size):
    x,y = [],[]
    x.append(xmin)
    y.append(h)
    x.append(xmax)
    y.append(h)
    x.append(xmax)
    y.append(h+size)
    x.append(xmin)
    y.append(h+size)
    x.append(xmin)
    y.append(h)
    return x,y

def line_chr(h,xmin,xmax,size):
    x,y = [],[]
    x.append(xmin)
    y.append(h)
    x.append(xmax)
    y.append(h)
    radius = size / 2
    # 计算两点的相对位置，确定半圆的起始角度和终止角度
    angle = np.arctan2(-size, 0)
    # 生成半圆的坐标点
    t = np.arange(angle, angle + np.pi, np.pi/180)
    x1 = radius * np.cos(t)
    y1 = radius * np.sin(t)
    x = x + list(xmax + x1)
    y = y + list(h+ radius + y1)
    x.append(xmax)
    y.append(h+size)
    x.append(xmin)
    y.append(h+size)
    # 计算两点的相对位置，确定半圆的起始角度和终止角度
    angle = np.arctan2(size, 0)
    # 生成半圆的坐标点
    t = np.arange(angle, angle + np.pi, np.pi/180)
    x1 = radius * np.cos(t)
    y1 = radius * np.sin(t)
    x = x + list(xmin + x1)
    y = y + list(h+ radius + y1)
    return x,y


def transform_pt(ch, pos, r, total_size,start_list):
    # convert chromosome position to axis coords
    rad = to_radian(pos + start_list[ch], total_size)
    return r*cos(rad), r*sin(rad)

def plot_colineartly_block(start_list,chrolist,data,radius, total_size,chr_reverse,lens):# 环图
    if chr_reverse == 'False':
        # pos1,pos2,pos3,pos4 = gff[id1]['end'],gff[id2]['end'],gff[id3]['end'],gff[id4]['end']
        pass
    else:
        if lens[chrolist[0]]['rev'] == '-':
            # pos1,pos2,pos3,pos4 = gff[id1]['end'],gff[id2]['end']
            data[0],data[2] = lens[chrolist[0]]['end'] - data[0],lens[chrolist[0]]['end'] - data[2]
        if lens[chrolist[1]]['rev'] == '-':
            # pos1,pos2,pos3,pos4 = gff[id1]['end'],gff[id2]['end']
            data[1],data[3] = lens[chrolist[1]]['end'] - data[0],lens[chrolist[1]]['end'] - data[2]

    # 1 -> 3 弧线
    t = np.arange(to_radian(min(data[0],data[2]) + start_list[chrolist[0]], total_size), to_radian(max(data[0],data[2]) + start_list[chrolist[0]], total_size), math.pi/720.)
    x, y = radius*np.cos(t), radius*np.sin(t)
    if data[0] < data[2]:
        x,y = list(x),list(y)
    else:
        x,y = list(x[::-1]),list(y[::-1])
    

    # 3 -> 4 贝塞尔曲线
    ex1x, ex1y = transform_pt(chrolist[0], data[2], radius, total_size,start_list)
    ex2x, ex2y = transform_pt(chrolist[1], data[3], radius, total_size,start_list)
    ratio = .5
    x0 = [ex1x, ex1x*ratio, ex2x*ratio, ex2x]
    y0 = [ex1y, ex1y*ratio, ex2y*ratio, ex2y]
    step = .01
    t = arange(0, 1+step, step)
    xt = Bezier(x0, t)
    yt = Bezier(y0, t)

    x = x + [i for i in xt]
    y = y + [i for i in yt]

    # 
    t = np.arange(to_radian(min(data[1],data[3]) + start_list[chrolist[1]], total_size), to_radian(max(data[1],data[3]) + start_list[chrolist[1]], total_size), math.pi/720.)
    x1, y1 = radius*np.cos(t), radius*np.sin(t)

    if data[1] > data[3]:
        x1,y1 = list(x1),list(y1)
    else:
        x1,y1 = list(x1[::-1]),list(y1[::-1])
    
    x = x + x1
    y = y + y1

    #
    ex1x, ex1y = transform_pt(chrolist[1], data[1], radius, total_size,start_list)
    ex2x, ex2y = transform_pt(chrolist[0], data[0], radius, total_size,start_list)
    ratio = .5
    x0 = [ex1x, ex1x*ratio, ex2x*ratio, ex2x]
    y0 = [ex1y, ex1y*ratio, ex2y*ratio, ex2y]
    step = .01
    t = arange(0, 1+step, step)
    xt = Bezier(x0, t)
    yt = Bezier(y0, t)

    x = x + [i for i in xt]
    y = y + [i for i in yt]

    x.append(x[0])
    y.append(y[0])
    return x,y



def readcolineartly1(file,blockmin,step,genepairsfile_type,gff,chr_list,genelist,p):# 块首尾 p=True 两个基因中有一个就算
    block = []
    focus_list = []
    if genepairsfile_type == 'MCScanX':
        alphagenepairs = open(file, 'r', encoding='utf-8').read()
        file = alphagenepairs.split('## Alignment ')[1:]
        for local in file:
            local = local.strip('\n').split('\n')[1:]
            if len(local) < blockmin:
                continue
            else:
                # print(local)
                data = [local[0].split('\t')[1],local[0].split('\t')[2],local[-1].split('\t')[1],local[-1].split('\t')[2]]
                if data[0] not in gff.keys() or data[1] not in gff.keys():
                    continue
                chr1,chr2 = gff[data[0]]['chr'],gff[data[1]]['chr']
                if chr1 not in chr_list or chr2 not in chr_list:
                    continue
                # print(local)
                for line in local:
                    lt = line.split('\t')
                    if p == 'True':
                        if lt[1] in genelist or lt[2] in genelist:
                            focus_list.append([lt[1],lt[2]])
                    else:
                        if lt[1] in genelist and lt[2] in genelist:
                            focus_list.append([lt[1],lt[2]])
                if chr1 == chr2 and min(abs(gff[data[0]]['order'] - gff[data[1]]['order']),abs(gff[data[0]]['order'] - gff[data[3]]['order']),abs(gff[data[2]]['order'] - gff[data[1]]['order']),abs(gff[data[2]]['order'] - gff[data[3]]['order'])) > step:
                    block.append(data)
                elif chr1 != chr2:
                    block.append(data)
    elif genepairsfile_type == 'WGDI':
        alphagenepairs = open(file, 'r', encoding='utf-8').read()
        file = alphagenepairs.split('# Alignment ')[1:]
        for local in file:
            local = local.strip('\n').split('\n')[1:]
            # print(local)
            if len(local) < blockmin:
                continue
            else:
                # print(local)
                data = [local[0].split()[0],local[0].split()[2],local[-1].split()[0],local[-1].split()[2]]
                if data[0] not in gff.keys() or data[1] not in gff.keys():
                    continue
                chr1,chr2 = gff[data[0]]['chr'],gff[data[1]]['chr']
                if chr1 not in chr_list or chr2 not in chr_list:
                    continue

                for line in local:
                    lt = line.split()
                    if p == 'True':
                        if lt[0] in genelist or lt[2] in genelist:
                            focus_list.append([lt[0],lt[2]])
                    else:
                        if lt[0] in genelist and lt[2] in genelist:
                            focus_list.append([lt[0],lt[2]])

                if chr1 == chr2 and min(abs(gff[data[0]]['order'] - gff[data[1]]['order']),abs(gff[data[0]]['order'] - gff[data[3]]['order']),abs(gff[data[2]]['order'] - gff[data[1]]['order']),abs(gff[data[2]]['order'] - gff[data[3]]['order'])) > step:
                    block.append(data)
                elif chr1 != chr2:
                    block.append(data)

    else:
        print('genepairsfile_type error: File Format not recognized! 共线性文件读取方式为block两头 2023.8.8修改')
        exit()
    return block,focus_list

def readcolineartly(file,blockmin,step,genepairsfile_type,gff,chr_list):# 块首尾
    block = []
    if genepairsfile_type == 'MCScanX':
        alphagenepairs = open(file, 'r', encoding='utf-8').read()
        file = alphagenepairs.split('## Alignment ')[1:]
        for local in file:
            local = local.strip('\n').split('\n')[1:]
            if len(local) < blockmin:
                continue
            else:
                # print(local)
                data = [local[0].split('\t')[1],local[0].split('\t')[2],local[-1].split('\t')[1],local[-1].split('\t')[2]]
                if data[0] not in gff.keys() or data[1] not in gff.keys():
                    continue
                chr1,chr2 = gff[data[0]]['chr'],gff[data[1]]['chr']
                # print(chr1,chr2,chr_list)

                if chr1 not in chr_list or chr2 not in chr_list:
                    # print('^^^^^^^^^^^^^^^^')
                    continue
                # print(local)
                if chr1 == chr2 and min(abs(gff[data[0]]['order'] - gff[data[1]]['order']),abs(gff[data[0]]['order'] - gff[data[3]]['order']),abs(gff[data[2]]['order'] - gff[data[1]]['order']),abs(gff[data[2]]['order'] - gff[data[3]]['order'])) > step:
                    block.append(data)
                elif chr1 != chr2:
                    block.append(data)
    elif genepairsfile_type == 'WGDI':
        alphagenepairs = open(file, 'r', encoding='utf-8').read()
        if "the\t# Alignment" in alphagenepairs:
            file = alphagenepairs.split('the\t')[1:]
        else:
            file = alphagenepairs.split('# Alignment ')[1:]
        # print(file[0])
        for local in file:
            local = local.strip('\n').split('\n')[1:]
            # print(local)
            if len(local) < blockmin:
                continue
            else:
                # print(local)
                data = [local[0].split()[0],local[0].split()[2],local[-1].split()[0],local[-1].split()[2]]
                # print(data)
                if data[0] not in gff.keys() or data[1] not in gff.keys():
                    continue
                chr1,chr2 = gff[data[0]]['chr'],gff[data[1]]['chr']
                # print(chr1,chr2)
                if chr1 not in chr_list or chr2 not in chr_list:
                    continue
                if chr1 == chr2 and min(abs(gff[data[0]]['order'] - gff[data[1]]['order']),abs(gff[data[0]]['order'] - gff[data[3]]['order']),abs(gff[data[2]]['order'] - gff[data[1]]['order']),abs(gff[data[2]]['order'] - gff[data[3]]['order'])) > step:
                    block.append(data)
                elif chr1 != chr2:
                    block.append(data)
    elif genepairsfile_type == 'ColinearScan':
        alphagenepairs = open(file, 'r', encoding='utf-8').read()
        file = alphagenepairs.split('# Alignment ')[1:]
        # print(file[0])
        for local in file:
            local = local.strip('\n').split('\n')[1:]
            # print(local)
            if len(local) < blockmin:
                continue
            else:
                # print(local)
                data = [local[0].split()[0],local[0].split()[2],local[-1].split()[0],local[-1].split()[2]]
                # print(data)
                if data[0] not in gff.keys() or data[1] not in gff.keys():
                    continue
                chr1,chr2 = gff[data[0]]['chr'],gff[data[1]]['chr']
                # print(chr1,chr2)
                if chr1 not in chr_list or chr2 not in chr_list:
                    continue
                if chr1 == chr2 and min(abs(gff[data[0]]['order'] - gff[data[1]]['order']),abs(gff[data[0]]['order'] - gff[data[3]]['order']),abs(gff[data[2]]['order'] - gff[data[1]]['order']),abs(gff[data[2]]['order'] - gff[data[3]]['order'])) > step:
                    block.append(data)
                elif chr1 != chr2:
                    block.append(data)
    else:
        print('genepairsfile_type error: File Format not recognized! 共线性文件读取方式为block两头 2023.8.8修改')
        exit()
    return block


def read_repeatgff3(gff3,chr_list):
    gff = pd.read_csv(gff3,header = None, sep='\t', comment='#')
    gff = gff.drop(gff[(gff[0].apply(lambda x: x not in chr_list))].index)# 条件过滤
    return gff

def read_fasta(file):
    fasta = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
    return fasta

def count_str(seq):
    str_dic = {'A':0,'T':0,'C':0,'G':0,'N':0}
    for i in seq:
        str_dic[i] = str_dic[i] + 1
    if str_dic['G'] + str_dic['C'] == 0:
        return 0
    else:
        return (str_dic['G'] + str_dic['C']) / len(seq)
def get_genomeGC(genome,windows,step,chr_list):
    GC = []
    for chro in genome:
        if chro not in chr_list:
            continue
        length = len(genome[chro].seq)
        fasta = str(genome[chro].seq).upper()
        job = int((length-windows)/step)
        for i in trange(job, desc=chro+' GC ->'):
            if i == job - 1:
                seq0 = fasta[i*step:]
                end = length
            seq0 = fasta[i*step:i*step+windows]
            end = i*step+windows
            GC_ = count_str(seq0)
            GC.append([chro,chro+'_'+str(i),i*step,end,'',i,'',GC_])
    data=pd.DataFrame(GC)
    return data


def count_Coverage(list0):
    overlap = list0.sum(axis=0)/len(list0)#平均每个位点上的重叠量
    list0 = np.int64(list0>0)
    Coverage = list0.sum(axis=0)/len(list0)#覆盖度
    return Coverage,overlap

def get_genedensity(gff,lens,windows,step,chr_list,a,b):
    Coverages,overlaps = [],[]
    gff = pd.read_csv(gff,header = None, sep='\t', comment='#', dtype=str)
    dic0 = gff.groupby(0).groups
    gff[a] = gff[a].astype(int)
    gff[b] = gff[b].astype(int)
    # print(gff)
    for key in dic0.keys():
        if key not in chr_list:
            continue
        local = gff.loc[dic0[key]].sort_values(by=[0,a],ascending= [True,True])
        local.reset_index(drop=True,inplace=True)
        chro = np.zeros(lens[key]['end'])
        for index,row in local.iterrows():
            chro[row[a]-1:row[b]] = chro[row[a]-1:row[b]] + 1
        length = lens[key]['end']
        job = int((length-windows)/step)
        for i in trange(job, desc=key+' density ->'):
            if i == job - 1:
                seq0 = chro[i*step:]
                end = length
            seq0 = chro[i*step:i*step+windows]
            end = i*step+windows
            Coverage,overlap = count_Coverage(seq0)
            Coverages.append([key,key+'_'+str(i),i*step,end,'',i,'',Coverage])
            overlaps.append([key,key+'_'+str(i),i*step,end,'',i,'',overlap])
    Coverages=pd.DataFrame(Coverages)
    overlaps=pd.DataFrame(overlaps)
    return Coverages,overlaps

def rad_to_coord(angle, radius):
    return radius*cos(angle), radius*sin(angle)

def plot_arc_local(start, stop, radius,color,alpha,size):
    # start, stop measured in radian
    t = arange(start, stop, pi/720.)
    x, y = radius*cos(t), radius*sin(t)
    plt.plot(x, y, color, alpha=alpha, linewidth=size)# 染色体圆弧

def drw_line(chr_list,s_e_dic,rl,text_b,rt,color,alpha,size):# 绘制线
    for i in chr_list:
        # print(i,'\n**********************************************')
        start = s_e_dic[i]['start']
        end = s_e_dic[i]['end']
        plot_arc_local(start, end, rl,color,alpha,size)
        if text_b:
            label_x, label_y = rad_to_coord((start+end)/2, rt)
            text(label_x, label_y, i, horizontalalignment="center", verticalalignment="center", fontsize = 7, color = 'black')
        else:
            pass

def readblast(genepairs,genepairsfile_type,blockl,class1):
    # print(block)
    blast_dic = {}
    one_gene = []
    alphagenepairs = open(genepairs, 'r', encoding='utf-8')
    if genepairsfile_type == 'famCircle':
        for row in alphagenepairs:
            if (row[0] == '\n'):
                continue
            elif (row[0] == '#'):
                lt = row.strip('\n').split(' ')
                block = {}
                for i in lt:
                    if '=' in str(i):
                        lt0 = i.split('=')
                        block[str(lt0[0])] = str(lt0[1])
                N = int(block['N'])
                if N >= int(blockl):
                    class1 = True
                else:
                    class1 = False
            else:
                if class1:
                    lt = row.strip('\n').split()
                    id1, id2 = str(lt[0]),str(lt[1]) 
                    if id1 not in one_gene:
                        one_gene.append(id1)
                        blast_dic[id1] = [id2]
                    else:
                        blast_dic[id1].append(id2)
    elif genepairsfile_type == 'WGDI':
        for row in alphagenepairs:
            if (row[0] == '\n'):
                continue
            elif (row[0] == '#'):
                lt = row.strip('\n').split(' ')
                block = {}
                for i in lt:
                    if '=' in str(i):
                        lt0 = i.split('=')
                        block[str(lt0[0])] = str(lt0[1])
                N = int(block['N'])
                # print(N)
                if N >= int(blockl):
                    class1 = True
                else :
                    class1 = False
            else:
                if class1:
                    lt = row.strip('\n').split()
                    id1, id2 = str(lt[0]),str(lt[2])
                    if id1 not in one_gene:
                        one_gene.append(id1)
                        blast_dic[id1] = [id2]
                    else:
                        blast_dic[id1].append(id2)
    elif genepairsfile_type == 'ColinearScan':
        for row in alphagenepairs:
            if (row[0] == '\n',row[0] == '+'):
                continue
            elif (row[:3] == 'the'):
                lt = row.strip('\n').split()
                N = int(lt[-1])
                if N >= int(blockl):
                    class1 = True
                else :
                    class1 = False
            else:
                if class1:
                    lt = row.strip('\n').split()
                    id1, id2 = str(lt[0]),str(lt[2])
                    if id1 not in one_gene:
                        one_gene.append(id1)
                        blast_dic[id1] = [id2]
                    else:
                        blast_dic[id1].append(id2)
    elif genepairsfile_type == 'MCScanX':
        for row in alphagenepairs:
            if (row[0] == '\n'):
                continue
            elif(row[:12] == '## Alignment'):
                lt = row.strip('\n').split(' ')
                block = {}
                for i in lt:
                    if '=' in str(i):
                        lt0 = i.split('=')
                        block[str(lt0[0])] = str(lt0[1])
                N = int(block['N'])
                # print(N)
                if N >= int(blockl):
                    class1 = True
                else :
                    class1 = False
            elif ('#' not in row):
                if class1:
                    lt = row.strip('\n').split()
                    # print(lt)
                    if len(lt) == 5:
                        id1, id2 = str(lt[2]),str(lt[3])
                    elif len(lt) == 4:
                        id1, id2 = str(lt[1]),str(lt[2])
                    elif len(lt) == 6:
                        id1, id2 = str(lt[3]),str(lt[4])
                    else:
                        print(row)
                        print('Parse error!')
                        exit()
                    if id1 not in one_gene:
                        one_gene.append(id1)
                        blast_dic[id1] = [id2]
                    else:
                        blast_dic[id1].append(id2)
    elif genepairsfile_type == 'BLAST':
        num = 0
        name0_list = []
        for row in alphagenepairs:
            if (row[0] == '\n'):
                continue
            elif ('#' not in row):
                lt = row.strip('\n').split()
                if lt[0] == lt[1]:
                    continue
                if str(lt[0]) not in name0_list:
                    name0_list.append(str(lt[0]))
                    num = 1
                else:
                    if num <= int(blockl):
                        pass
                    else:
                        continue
                id1, id2 = str(lt[0]),str(lt[1])
                if id1 not in one_gene:
                    one_gene.append(id1)
                    blast_dic[id1] = [id2]
                else:
                    blast_dic[id1].append(id2)
    else:
        print('genepairsfile_type error: File Format not recognized!')
        exit()
    return blast_dic

def readcolineartly_dotplot(file,blockmin,step,genepairsfile_type,gff,chr_list):# 块首尾
    block = []
    if genepairsfile_type == 'MCScanX':
        alphagenepairs = open(file, 'r', encoding='utf-8').read()
        file = alphagenepairs.split('## Alignment ')[1:]
        for local in file:
            local = local.strip('\n').split('\n')[1:]
            if len(local) < blockmin:
                continue
            else:
                data = [[i.split('\t')[1],i.split('\t')[2]] for i in local]
                if data[0][0] not in gff.keys() or data[0][1] not in gff.keys():
                    continue
                chr1,chr2 = gff[data[0][0]]['chr'],gff[data[0][1]]['chr']
                if chr1 not in chr_list or chr2 not in chr_list:
                    continue
                if chr1 == chr2 and min(abs(gff[data[0][0]]['order'] - gff[data[0][1]]['order']),abs(gff[data[0][0]]['order'] - gff[data[-1][1]]['order']),abs(gff[data[-1][0]]['order'] - gff[data[0][1]]['order']),abs(gff[data[-1][0]]['order'] - gff[data[-1][1]]['order'])) > step:
                    block.append(data)
                elif chr1 != chr2:
                    block.append(data)
    elif genepairsfile_type == 'WGDI':
        alphagenepairs = open(file, 'r', encoding='utf-8').read()
        file = alphagenepairs.split('# Alignment ')[1:]
        for local in file:
            local = local.strip('\n').split('\n')[1:]
            if len(local) < blockmin:
                continue
            else:
                data = [[i.split()[0],i.split()[2]] for i in local]
                if data[0][0] not in gff.keys() or data[0][1] not in gff.keys():
                    continue
                chr1,chr2 = gff[data[0][0]]['chr'],gff[data[0][1]]['chr']
                if chr1 not in chr_list or chr2 not in chr_list:
                    continue
                if chr1 == chr2 and min(abs(gff[data[0][0]]['order'] - gff[data[0][1]]['order']),abs(gff[data[0][0]]['order'] - gff[data[-1][1]]['order']),abs(gff[data[-1][0]]['order'] - gff[data[0][1]]['order']),abs(gff[data[-1][0]]['order'] - gff[data[-1][1]]['order'])) > step:
                    block.append(data)
                elif chr1 != chr2:
                    block.append(data)

    else:
        print('genepairsfile_type error: File Format not recognized! 共线性文件读取方式为完整block 2023.8.21修改')
        exit()
    return block

def cat_gene(gene1,gene2):
    return '_'.join([gene1,gene2])

import time
from tqdm import tqdm
def readblast_large(genepairs_,blast_reverse,score_levels_,position_,levels_,tandem_,gff,lens,chr1_start,chr2_start,step1,step2):
    print("Count the number of lines in the BLAST file")
    start_time = time.time()  # 记录开始时间
    total_lines = sum(1 for _ in open(genepairs_,'r'))
    end_time = time.time()  # 记录结束时间
    print(f"Function execution time: {end_time - start_time:.2f} seconds")
    x1,x2,x3,y1,y2,y3 = [],[],[],[],[],[]
    # 初始化进度条

    with tqdm(total=math.ceil(total_lines/100000), desc="Reading BLAST", unit="block") as pbar:
        # 分块读取文件
        chunk_iter = pd.read_csv(genepairs_,header = None, sep='\t', comment='#', chunksize=100000)
        processed_lines = 0
        tail = pd.DataFrame()
        for chunk in chunk_iter:
            processed_lines += len(chunk)
            chunk.loc[:, 12] = chunk.apply(lambda x : cat_gene(x[0],x[1]), axis=1)
            chunk = chunk.drop_duplicates(subset=12 , keep='first')
            del chunk[12]
            if blast_reverse == 'True':
                chunk[[0, 1]] = chunk[[1, 0]]
            chunk = chunk.drop(chunk[(chunk[0] == chunk[1])].index)#.sort_values(by=[0,11],ascending= [True,False])
            last_key = list(chunk[0].to_list())[-1]
            if not tail.empty:
                chunk = pd.concat([tail, chunk], ignore_index=True)
                chunk.reset_index(drop=True,inplace=True)
            dic = chunk.groupby(0).groups# 按照第几列分组
            # 判断是否为最后一块
            if processed_lines == total_lines:
                # print("This is the last chunk.")
                tail = pd.DataFrame()
            else:
                tail = chunk.loc[dic[last_key]]
                tail.reset_index(drop=True,inplace=True)
                del dic[last_key]

            for key in dic.keys():
                # print(" ******** ",key," ******** ")
                local = chunk.loc[dic[key]]
                local.reset_index(drop=True,inplace=True)

                if key not in gff.keys():
                    continue
                if score_levels_ == 'False':
                    if tandem_ == 'True':
                        object0 = local[1].to_list()
                        object0 = list(dict.fromkeys(object0))# 去除tandem保留顺序
                        for i in range(levels_[0]):
                            if i < len(object0):
                                if object0[i] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[object0[i]]['chr'] not in lens.keys():
                                    continue
                                if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                    x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                                else:
                                    x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                                if 'rev' in lens[gff[object0[i]]['chr']].keys() and lens[gff[object0[i]]['chr']]['rev'] == '-':
                                    y = -(chr2_start[gff[object0[i]]['chr']] + lens[gff[object0[i]]['chr']][position_] - gff[object0[i]][position_] + 1)*step2
                                else:
                                    y = -(chr2_start[gff[object0[i]]['chr']] + gff[object0[i]][position_])*step2
                                x1.append(x)
                                y1.append(y)
                        for i in range(levels_[0],levels_[0]+levels_[1]):
                            if i < len(object0):
                                if object0[i] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[object0[i]]['chr'] not in lens.keys():
                                    continue
                                if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                    x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                                else:
                                    x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                                if 'rev' in lens[gff[object0[i]]['chr']].keys() and lens[gff[object0[i]]['chr']]['rev'] == '-':
                                    y = -(chr2_start[gff[object0[i]]['chr']] + lens[gff[object0[i]]['chr']][position_] - gff[object0[i]][position_] + 1)*step2
                                else:
                                    y = -(chr2_start[gff[object0[i]]['chr']] + gff[object0[i]][position_])*step2
                                x2.append(x)
                                y2.append(y)
                        for i in range(levels_[0]+levels_[1],levels_[0]+levels_[1]+levels_[2]):
                            if i < len(object0):
                                if object0[i] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[object0[i]]['chr'] not in lens.keys():
                                    continue
                                if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                    x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                                else:
                                    x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                                if 'rev' in lens[gff[object0[i]]['chr']].keys() and lens[gff[object0[i]]['chr']]['rev'] == '-':
                                    y = -(chr2_start[gff[object0[i]]['chr']] + lens[gff[object0[i]]['chr']][position_] - gff[object0[i]][position_] + 1)*step2
                                else:
                                    y = -(chr2_start[gff[object0[i]]['chr']] + gff[object0[i]][position_])*step2
                                x3.append(x)
                                y3.append(y)
                    else:
                        object0 = local[1].to_list()
                        for i in range(levels_[0]):
                            if i < len(object0):
                                if object0[i] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[object0[i]]['chr'] not in lens.keys():
                                    continue
                                if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                    x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                                else:
                                    x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                                if 'rev' in lens[gff[object0[i]]['chr']].keys() and lens[gff[object0[i]]['chr']]['rev'] == '-':
                                    y = -(chr2_start[gff[object0[i]]['chr']] + lens[gff[object0[i]]['chr']][position_] - gff[object0[i]][position_] + 1)*step2
                                else:
                                    y = -(chr2_start[gff[object0[i]]['chr']] + gff[object0[i]][position_])*step2
                                x1.append(x)
                                y1.append(y)
                        for i in range(levels_[0],levels_[0]+levels_[1]):
                            if i < len(object0):
                                if object0[i] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[object0[i]]['chr'] not in lens.keys():
                                    continue
                                if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                    x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                                else:
                                    x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                                if 'rev' in lens[gff[object0[i]]['chr']].keys() and lens[gff[object0[i]]['chr']]['rev'] == '-':
                                    y = -(chr2_start[gff[object0[i]]['chr']] + lens[gff[object0[i]]['chr']][position_] - gff[object0[i]][position_] + 1)*step2
                                else:
                                    y = -(chr2_start[gff[object0[i]]['chr']] + gff[object0[i]][position_])*step2
                                x2.append(x)
                                y2.append(y)
                        for i in range(levels_[0]+levels_[1],levels_[0]+levels_[1]+levels_[2]):
                            if i < len(object0):
                                if object0[i] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[object0[i]]['chr'] not in lens.keys():
                                    continue
                                if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                    x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                                else:
                                    x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                                if 'rev' in lens[gff[object0[i]]['chr']].keys() and lens[gff[object0[i]]['chr']]['rev'] == '-':
                                    y = -(chr2_start[gff[object0[i]]['chr']] + lens[gff[object0[i]]['chr']][position_] - gff[object0[i]][position_] + 1)*step2
                                else:
                                    y = -(chr2_start[gff[object0[i]]['chr']] + gff[object0[i]][position_])*step2
                                x3.append(x)
                                y3.append(y)
                else:
                    score_levels = [float(i) for i in score_levels_.split(',')]
                    for index,row in local.iterrows():
                        if row[2] >= score_levels[0]:
                            if row[1] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[row[1]]['chr'] not in lens.keys():
                                continue
                            if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                            else:
                                x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                            if 'rev' in lens[gff[row[1]]['chr']].keys() and lens[gff[row[1]]['chr']]['rev'] == '-':
                                y = -(chr2_start[gff[row[1]]['chr']] + lens[gff[row[1]]['chr']][position_] - gff[row[1]][position_] + 1)*step2
                            else:
                                y = -(chr2_start[gff[row[1]]['chr']] + gff[row[1]][position_])*step2
                            x1.append(x)
                            y1.append(y)
                        elif row[2] >= score_levels[1]:
                            if row[1] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[row[1]]['chr'] not in lens.keys():
                                continue
                            if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                            else:
                                x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                            if 'rev' in lens[gff[row[1]]['chr']].keys() and lens[gff[row[1]]['chr']]['rev'] == '-':
                                y = -(chr2_start[gff[row[1]]['chr']] + lens[gff[row[1]]['chr']][position_] - gff[row[1]][position_] + 1)*step2
                            else:
                                y = -(chr2_start[gff[row[1]]['chr']] + gff[row[1]][position_])*step2
                            x2.append(x)
                            y2.append(y)
                        elif row[2] >= score_levels[2]:
                            if row[1] not in gff.keys() or gff[key]['chr'] not in lens.keys() or gff[row[1]]['chr'] not in lens.keys():
                                continue
                            if 'rev' in lens[gff[key]['chr']].keys() and lens[gff[key]['chr']]['rev'] == '-':
                                x = (chr1_start[gff[key]['chr']] + lens[gff[key]['chr']][position_] - gff[key][position_] + 1)*step1
                            else:
                                x = (chr1_start[gff[key]['chr']] + gff[key][position_])*step1
                            if 'rev' in lens[gff[row[1]]['chr']].keys() and lens[gff[row[1]]['chr']]['rev'] == '-':
                                y = -(chr2_start[gff[row[1]]['chr']] + lens[gff[row[1]]['chr']][position_] - gff[row[1]][position_] + 1)*step2
                            else:
                                y = -(chr2_start[gff[row[1]]['chr']] + gff[row[1]][position_])*step2
                            x3.append(x)
                            y3.append(y)
            # pbar.update(len(chunk))
            pbar.update()
            
    return x1,x2,x3,y1,y2,y3


def value_to_color(value, cmap, vmin, vmax):
    norm = Normalize(vmin=vmin, vmax=vmax)
    return cmap(norm(value))

import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
def make_cmap(seabornc):
    # 创建渐变颜色映射
    colors = sns.color_palette(seabornc, as_cmap=True)
    cmap = LinearSegmentedColormap.from_list("Custom", colors(np.linspace(0, 1, 256)))
    return cmap


def write_fasta(seq0,id0,file):# 写入fasta格式
    rec = SeqRecord(Seq((str(seq0))), id=id0)
    file = open(file, "a+")# 切记追加写入
    SeqIO.write(rec,file,"fasta")
    file.close()

def align(align_software,path,inputf,outputf):
    # if align_software == 'mafft':
    #     mafft_cline = MafftCommandline(
    #         cmd=path, input=inputf, auto=True)
    #     stdout, stderr = mafft_cline()
    #     align = AlignIO.read(StringIO(stdout), "fasta")
    #     AlignIO.write(align, outputf, "fasta")
    # if align_software == 'muscle':
    #     muscle_cline = MuscleCommandline(
    #         cmd=path, input=inputf, out=outputf, seqtype="protein", clwstrict=True)
    #     stdout, stderr = muscle_cline()
    # if align_software == "clustalw":
    #     Clustalw_cline = ClustalwCommandline(cmd = path, infile=inputf, outfile=outputf)
    #     stdout, stderr = Clustalw_cline()

    if align_software == 'mafft':
        command = [path, '--auto', inputf]
        result = subprocess.run(command, capture_output=True, text=True)
        align = AlignIO.read(StringIO(result.stdout), "fasta")
        AlignIO.write(align, outputf, "fasta")
    elif align_software == 'muscle':
        command = [path, '-in', inputf, '-out', outputf]
        result = subprocess.run(command, capture_output=True, text=True)
    elif align_software == 'clustalw':
        command = [path, '-INFILE='+inputf, '-OUTPUT=FASTA', '-OUTFILE='+outputf]
        result = subprocess.run(command, capture_output=True, text=True)
    else:
        print('多序列比对软件支持muscle/mafft/clustalw 请重新选择！')
        exit()

# , align=True, outorder="ALIGNED", convert=True, output="pir", clwstrictout=outputf


# def plot_cap(angle, clockwise):
#     radius=sm_radius
#     # angle measured in radian, clockwise is boolean 鸭舌
#     if clockwise: 
#         t = np.arange(angle, angle+math.pi, math.pi/30.)
#     else: 
#         t = np.arange(angle, angle-math.pi, -math.pi/30.)
#     x, y = radius*np.cos(t), radius*np.sin(t)
#     middle_r = (radius_a+radius_b)/2
#     x, y = x + middle_r*np.cos(angle), y + middle_r*np.sin(angle)
#     plot(x, y, "k-", alpha=.5)



###################################################################

# def read_ks1(file0):
#     blast0=pd.read_csv(file0,header = 0, sep='\t')
#     blast0["all"] = blast0["id1"].map(str)+"%"+blast0["id2"].map(str)
#     blast0["all_"] = blast0["all"].apply(scale)
#     df = blast0.drop_duplicates("all_",keep='first')
#     ks_ = df[["all","ks_NG86"]]# 保留两列
#     del blast0,df
#     gc.collect()
#     ks = ks_.set_index("all").T.to_dict()
#     lt = list(ks.keys())
#     for i in trange(len(lt)):
#         key = lt[i]
#         ks[key] = float(ks[key]["ks_NG86"])
#     return ks

# def read_family1(file):
#     family = []
#     family_dic = {}
#     dic_motif = {}
#     for i in open(file, 'r', encoding='utf-8'):
#         if i[0] != '#' or i[0] != '\n':
#             line = i.strip('\n').split()
#             if line[0] not in family:
#                 family.append(line[0])
#                 if len(line) > 1:
#                     family_dic[line[0]] = line[1]
#                     if line[2] not in dic_motif.keys():
#                         dic_motif[line[2]] = line[1]
#                 else:
#                     family_dic[line[0]] = '#4d5aaf'
#     return family,family_dic,dic_motif
# def read_gff(gff):
#     gene = {}
#     gff = open(gff, 'r', encoding='utf-8').readlines()
#     for j in range(len(gff)):
#         i = gff[j]
#         if i[0] != '#' or len(i.strip('\n')) != '':
#             line = i.strip('\n').split()
#             dic0 = {}
#             dic0['chr'] = line[0]
#             dic0['end'] = int(line[2])
#             dic0['order'] = int(line[5])
#             gene[line[1]] = dic0
#     return gene

# def read_lens(lens):
#     chr_dic = {}
#     chr_list = []
#     for i in open(lens, 'r', encoding='utf-8'):
#         if i[0] != '#' or len(i.strip('\n')) != '':
#             line = i.strip('\n').split()
#             dic0 = {}
#             dic0['end'] = int(line[1])
#             dic0['order'] = int(line[2])
#             chr_dic[line[0]] = dic0
#             chr_list.append(line[0])
#     return chr_dic,chr_list



def qiepian0(name):
    return str(name).split("%")[0]
def qiepian1(name):
    return str(name).split("%")[1]



def read_blast0(blast0):# 矩阵方法
    print("读取BLAST")
    df1=pd.read_csv(blast0,header = None, sep='\t')
    df = df1[[0,1]]# 保留两列
    # print(df)
    df.columns=["Query","Subject"]# 添加列名
    # print(df.groupby("Query").size())# 分组信息
    df = df.drop(df[df["Query"] == df["Subject"]].index)
    data_dict = df.groupby("Query").Subject.apply(list).to_dict()# 分组后转列表
    return data_dict

def scale(name):
    lt = name.split("%")
    lt.sort()
    name = "".join(lt)
    return name

def read_ks0(file0):
    ks = {}
    file = open(file0, 'r', encoding='utf-8').readlines()
    for j in trange(len(file), desc='read ks ->'):
        i = file[j]
        if i[0] == '#' or i[0] == '\n':
            continue
        elif 'id1' in i:
            continue
        else:
            line = i.strip('\n').split()
            if len(line) <= 4:
                continue
            ks[str(line[0]) + '_' + str(line[1])] = float(line[3])
    return ks

def read_WGDI0(file0,block):
    blast = {}
    class1 = False
    file = open(file0, 'r', encoding='utf-8').readlines()
    for j in trange(len(file), desc='read WGDI ->'):
        i = file[j]
        if i[0] == '\n':
            continue
        elif i[0] == '#':
            lt = i.strip('\n').split()
            for x in lt:
                if 'N=' in x:
                    length = x[2:]
            if int(length) < int(block):
                class1 = False
            else:
                class1 = True
        else:
            if class1:
                line = i.strip('\n').split()
                if line[0] not in blast.keys():
                    lt = []
                    lt.append(line[2])
                    blast[line[0]] = lt
                else:
                    blast[line[0]].append(line[2])
    return blast

def read_family(file):
    family = []
    for i in open(file, 'r', encoding='utf-8'):
        if i[0] != '#' or i[0] != '\n':
            line = i.strip('\n').split()
            if line[0] not in family:
                family.append(line[0])
    return family



def getRGB(dWave,maxPix=1,gamma=1):
    #dWave为波长；maxPix为最大值；gamma为调教参数
    waveArea = [380,440,490,510,580,645,780]
    minusWave = [0,440,440,510,510,645,780]
    deltWave = [1,60,50,20,70,65,35]
    for p in range(len(waveArea)):
        if dWave<waveArea[p]:
            break
    pVar = abs(minusWave[p]-dWave)/deltWave[p]
    rgbs = [[0,0,0],[pVar,0,1],[0,pVar,1],[0,1,pVar],
            [pVar,1,0],[1,pVar,0],[1,0,0],[0,0,0]]
    #在光谱边缘处颜色变暗
    if (dWave>=380) & (dWave<420):
        alpha = 0.3+0.7*(dWave-380)/(420-380)
    elif (dWave>=420) & (dWave<701):
        alpha = 1.0
    elif (dWave>=701) & (dWave<780):
        alpha = 0.3+0.7*(780-dWave)/(780-700)
    elif dWave==780:
        alpha = 0.3+0.7*(0.001)/(780-700)
    else:
        alpha = 0       #非可见区
    return [maxPix*(c*alpha)**gamma for c in rgbs[p]]

def return_col1(model,node,ks,min_,max_):
    # ks越小波长越长颜色越偏红，紫短红长
    if model == 'gradual':
        if ks <= 10:
            wavelength = 760-(((ks-min_)/(max_-min_))*(760-400))
        else:
            wavelength = 400+(((ks-min_)/(max_-min_))*(760-400))
        return getRGB(wavelength)
    elif model == 'jump':
        # col_l = ['#d7003a','#2ca9e1','#b8d200','#f6ad49']
        col_l = ["#A90C38","#B4193A","#BF233C","#CB2C3E","#D63440","#E13E43","#E64E4B","#EB5C53","#F06A5B","#F57763","#F88570","#F59885","#F1A99A","#ECBBB0","#E4CCC6","#D6D0D1","#C3C8D2","#B0C0D2","#9CB8D3","#86B1D3","#78A8CE","#6F9FC6","#6696BF","#5E8DB8","#5584B0","#4D7CA8","#4573A0","#3E6B97","#36628F","#2E5A87"]
        lt = node.split(',')
        col_l = [col_l[i*int(len(col_l)/(len(lt)-1))] for i in range(len(lt)-1)]
        # print(col_l)
        if ks <= 10:
            for i in range((len(lt)-1)):
                if ks < float(lt[i]):
                    return col_l[i-1]
            return col_l[(len(lt)-1)-1]

        else:
            for i in range((len(lt)-1)):
                if ks >= float(lt[i]):
                    return col_l[-i]
            return col_l[(len(lt)-1)+1]

def return_col(ks,min_,max_):
    # ks越小波长越长颜色越偏红，紫短红长
    if ks <= 10:
        wavelength = 760-(((ks-min_)/(max_-min_))*(760-400))
    else:
        wavelength = 400+(((ks-min_)/(max_-min_))*(760-400))
    return getRGB(wavelength)

def drawSpec1(min_,max_,gs):#min_,max_
    if max_ <= 10:
        ax1 = plt.subplot(gs)
        #波长越大ks越小
        pic = np.zeros([10,360,3])
        rgb = [getRGB(d) for d in range(400,760)]
        pic = pic+rgb
        pic0 = np.zeros([360,10,3])
        for i in range(10):
            for j in range(360):
                pic0[j][i] = pic[i][j]
        plt.imshow(pic)
        plt.xticks(range(0,360,50),[str(format(max_ - i*((max_-min_)/8), '.4f')) for i in range(8)], rotation=30, fontsize=fontsize)      #y坐标轴ks越小颜色越红
        plt.yticks([])
    else:
        ax1 = plt.subplot(gs)
        #波长越大ks越小
        pic = np.zeros([10,360,3])
        rgb = [getRGB(d) for d in range(760,400,-1)]
        pic = pic+rgb
        pic0 = np.zeros([360,10,3])
        for i in range(10):
            for j in range(360):
                pic0[j][i] = pic[i][j]
        plt.imshow(pic)
        plt.xticks(range(0,360,50),[str(format(max_ - i*((max_-min_)/8), '.4f')) for i in range(8)], rotation=30, fontsize=fontsize)      #y坐标轴ks越小颜色越红
        plt.yticks([])

def drawSpec11(node,gs,fontsize):#
    col_l = ["#A90C38","#B4193A","#BF233C","#CB2C3E","#D63440","#E13E43","#E64E4B","#EB5C53","#F06A5B","#F57763","#F88570","#F59885","#F1A99A","#ECBBB0","#E4CCC6","#D6D0D1","#C3C8D2","#B0C0D2","#9CB8D3","#86B1D3","#78A8CE","#6F9FC6","#6696BF","#5E8DB8","#5584B0","#4D7CA8","#4573A0","#3E6B97","#36628F","#2E5A87"]
    # col_l = ['#d7003a','#ee7800','#fcc800','#3eb370','#165e83','#00a497','#bc64a4','#0d0015']
    lt = node.split(',')
    col_l = [col_l[i*int(len(col_l)/(len(lt)-1))] for i in range(len(lt)-1)]
    ax1 = plt.subplot(gs)
    for i in range(len(lt)):
        c = col_l[i]
        rect=plt.Rectangle(
                (0.1+(i*(0.8/(len(lt)-1))), 0.1),  # (x,y)矩形左下角
                0.8/(len(lt)-1),  # width长
                0.1,  # height宽
                color=c)
        ax1.add_patch(rect)
        plt.text(0.1+(i*(0.8/(len(lt)-1))), 0.1, str(lt[i]),fontsize=fontsize, ha='left',va='center')
    plt.axis('off')

def drawSpec_1(node,gs,fontsize,font_prop):#
    col_l = ["#A90C38","#B4193A","#BF233C","#CB2C3E","#D63440","#E13E43","#E64E4B","#EB5C53","#F06A5B","#F57763","#F88570","#F59885","#F1A99A","#ECBBB0","#E4CCC6","#D6D0D1","#C3C8D2","#B0C0D2","#9CB8D3","#86B1D3","#78A8CE","#6F9FC6","#6696BF","#5E8DB8","#5584B0","#4D7CA8","#4573A0","#3E6B97","#36628F","#2E5A87"]
    # col_l = ['#d7003a','#ee7800','#fcc800','#3eb370','#165e83','#00a497','#bc64a4','#0d0015']
    lt = node.split(',')
    col_l = [col_l[i*int(len(col_l)/(len(lt)-1))] for i in range(len(lt)-1)]
    # print(col_l)
    ax1 = plt.subplot(gs)
    for i in range(len(lt)-1):
        c = col_l[i]
        rect=plt.Rectangle(
                (0.1,(i*(0.8/(len(lt)-1)))+0.1),  # (x,y)矩形左下角
                0.1,  # width长
                0.8/(len(lt)-1),  # height宽
                color=c)
        ax1.add_patch(rect)
        plt.text(0.3,(i*(0.8/(len(lt)-1)))+0.1, str(lt[i]), ha='left',va='center', fontproperties=font_prop,fontsize=fontsize)
    plt.text(0.3,((len(lt)-1)*(0.8/(len(lt)-1)))+0.1, str(lt[-1]), ha='left',va='center', fontproperties=font_prop,fontsize=fontsize)
    plt.axis('off')

def drawSpec(min_,max_,gs,fontsize,font_prop):#min_,max_
    if max_ <= 10:
        ax1 = plt.subplot(gs)
        #波长越大ks越小
        pic = np.zeros([10,360,3])
        rgb = [getRGB(d) for d in range(400,760)]
        pic = pic+rgb
        pic0 = np.zeros([360,10,3])
        for i in range(10):
            for j in range(360):
                pic0[j][i] = pic[i][j]
        plt.imshow(pic0)
        plt.yticks(range(0,360,50),[str(format(max_ - i*((max_-min_)/8), '.4f')) for i in range(8)], rotation=30, fontproperties=font_prop, fontsize=fontsize)      #y坐标轴ks越小颜色越红
        plt.xticks([])
    else:
        ax1 = plt.subplot(gs)
        #波长越大ks越小
        pic = np.zeros([10,360,3])
        rgb = [getRGB(d) for d in range(760,400,-1)]
        pic = pic+rgb
        pic0 = np.zeros([360,10,3])
        for i in range(10):
            for j in range(360):
                pic0[j][i] = pic[i][j]
        plt.imshow(pic0)
        plt.yticks(range(0,360,50),[str(format(max_ - i*((max_-min_)/8), '.4f')) for i in range(8)], rotation=30, fontproperties=font_prop, fontsize=fontsize)      #y坐标轴ks越小颜色越红
        plt.xticks([])

def readname(name,chrolist):
    name = name.lower()
    chrolist = chrolist.lower()
    if ('^' in name):
        return name
    elif ('g' in name):
        if 'g' not in chrolist:
            name = name.replace('g', '^')
        else :
            lt = name.split('g')
            name = ''
            for i in range(len(lt)):
                if i < len(lt) - 2:
                    name = name + str(lt[i]) + 'g'
                elif i == len(lt) - 2:
                    name = name + str(lt[-2])
                else:
                    name = name + "^" + str(lt[-1])
        return name

# def config():
#     conf = configparser.ConfigParser()
#     conf.read(os.path.join(famCircle.__path__[0], 'conf.ini'))
#     return conf.items('ini')

# def load_conf(file, section):
#     conf = configparser.ConfigParser()
#     conf.read(file)
#     return conf.items(section)

def config():
    conf = configparser.ConfigParser()
    conf.read(os.path.join(famCircle.__path__[0], 'conf.ini'), encoding='utf-8')
    return conf.items('ini')

def load_conf(file, section):
    conf = configparser.ConfigParser()
    conf.read(file, encoding='utf-8')
    return conf.items(section)

def read_coliearity0(file,chrolist):
    f = open(file,'r', encoding='utf-8')
    pair_list = []
    for line in f:
        if line[0] == '#':
            continue
        else:
            lt = line.strip('\n').split()
            genepair = ''
            gene1 = readname(str(lt[0]),chrolist)
            gene2 = readname(str(lt[2]),chrolist)
            if gene1 < gene2:
               genepair = gene1 + " " + gene2
            else:
               genepair = gene2 + " " + gene1
            pair_list.append(genepair)
    return pair_list

def read_blast(file,chrolist):
    f = open(file,'r', encoding='utf-8')
    pair_list = []
    for line in f:
        if line[0] == '#':
            continue
        else:
            lt = line.strip('\n').split()
            genepair = ''
            gene1 = readname(str(lt[0]),chrolist)
            gene2 = readname(str(lt[1]),chrolist)
            if gene1 < gene2:
               genepair = gene1 + " " + gene2
            else:
               genepair = gene2 + " " + gene1
            pair_list.append(genepair)
    return pair_list

def calculate_coef(p0, p1, p2, p3):
    c = 3*(p1 - p0)
    b = 3*(p2 - p1) -c
    a = p3 - p0 - c - b
    return c, b, a

def Bezier(plist, t):
    # p0 : origin, p1, p2 :control, p3: destination
    p0, p1, p2, p3 = plist
    # calculates the coefficient values
    c, b, a = calculate_coef(p0, p1, p2, p3)
    tsquared = t**2
    tcubic = tsquared*t
    return a*tcubic + b*tsquared + c*t + p0

def plot_bezier_curve(ex1x, ex1y, ex2x, ex2y):
    ratio = .2
    dx,dy = ex2x-ex1x,ex2y-ex1y
    x = [ex1x, ex1x+(dx*ratio), ex2x-(dx*ratio), ex2x]
    y = [ex1y, ex1y+(dy*(1-ratio)), ex2y-(dy*(1-ratio)), ex2y]
    step = .01
    t = arange(0, 1+step, step)
    xt = Bezier(x, t)# 贝塞尔曲线
    yt = Bezier(y, t)
    return xt,yt

def gene_length(gfffile):
    # 读取基因长度，得到平均长度和最小最大长度
    f = open(gfffile,'r', encoding='utf-8')
    genelength = {}
    for row in f:
        if row[0] != '\n' and row[0] != '#':
            row = row.strip('\n').split('\t')
            if str(row[1]) in genelength.keys():
                continue
            if 'chr' == str(row[1])[3:6]:
                continue
            length = abs(int(row[3]) - int(row[2]))
            genelength[str(row[1])] = length
    f.close()
    lt = []
    for i in genelength.values():
        lt.append(i)
    pj = sum(lt)/len(lt)
    return pj

def gene_length_order(lensfile):
    # 读取基因长度，得到平均长度和最小最大长度
    f = open(lensfile,'r', encoding='utf-8')
    gene_sum = 0
    gene_index = 0
    for row in f:
        if row[0] != '\n' and row[0] != '#':
            row = row.strip('\n').split('\t')
            gene_sum = gene_sum + int(row[2])
            gene_index = gene_index + int(row[1])
    f.close()
    pj = gene_index/gene_sum
    return pj

def cds_to_pep(cds_file, pep_file, fmt='fasta'):
    records = list(SeqIO.parse(cds_file, fmt))
    for k in records:
        k.seq = k.seq.translate()
    SeqIO.write(records, pep_file, 'fasta')
    return True

def read_colinearscan(file):
    data, b, flag, num = [], [], 0, 1
    with open(file) as f:
        for line in f.readlines():
            line = line.strip()
            if re.match(r"the", line):
                num = re.search(r'\d+', line).group()
                b = []
                flag = 1
                continue
            if re.match(r"\>LOCALE", line):
                flag = 0
                p = re.split(':', line)
                if len(b) > 0:
                    data.append([num, b, p[1]])
                b = []
                continue
            if flag == 1:
                a = re.split(r"\s", line)
                b.append(a)
    return data


def read_mcscanx(fn):
    f1 = open(fn)
    data, b = [], []
    flag, num = 0, 0
    for line in f1.readlines():
        line = line.strip()
        if re.match(r"## Alignment", line):
            flag = 1
            if len(b) == 0:
                arr = re.findall(r"[\d+\.]+", line)[0]
                continue
            data.append([num, b, 0])
            b = []
            num = re.findall(r"\d+", line)[0]
            continue
        if flag == 0:
            continue
        a = re.split(r"\:", line)
        c = re.split(r"\s+", a[1])
        b.append([c[1], c[1], c[2], c[2]])
    data.append([num, b, 0])
    return data


def read_jcvi(fn):
    f1 = open(fn)
    data, b = [], []
    num = 1
    for line in f1.readlines():
        line = line.strip()
        if re.match(r"###", line):
            if len(b) == 0:
                continue
            data.append([num, b, 0])
            b = []
            num += 1
            continue
        a = re.split(r"\t", line)
        b.append([a[0], a[0], a[1], a[1]])
    data.append([num, b, 0])
    return data


def read_coliearity(fn):
    f1 = open(fn)
    data, b = [], []
    flag, num = 0, 0
    for line in f1.readlines():
        line = line.strip()
        if re.match(r"# Alignment", line):
            flag = 1
            if len(b) == 0:
                arr = re.findall(r'[\.\d+]+', line)
                continue
            data.append([arr[0], b, arr[2]])
            b = []
            arr = re.findall(r'[\.\d+]+', line)
            continue
        if flag == 0:
            continue
        b.append(re.split(r"\s", line))
    data.append([arr[0], b, arr[2]])
    return data

def read_ks(file, col):
    ks = pd.read_csv(file, sep='\t')
    ks.drop_duplicates(subset=['id1', 'id2'], keep='first', inplace=True)
    ks[col] = ks[col].astype(float)
    ks = ks[ks[col] >= 0]
    ks.index = ks['id1']+','+ks['id2']
    return ks[col]
