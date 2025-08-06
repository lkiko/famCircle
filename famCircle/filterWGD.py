# -*- encoding: utf-8 -*-
'''
@File        :filterWGD.py
@Time        :2021/09/28 11:21:18
@Author        :charles kiko
@Version        :1.0
@Contact        :charles_kiko@163.com
@Desc        :全基因的ks
'''


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from famCircle.bez import *
import sys
from scipy.stats.kde import gaussian_kde
from collections import OrderedDict
plt.rc('font',family='Arial')


class filterWGD():
    def __init__(self, options):
        self.area1 = "0,1000"
        self.area2 = "0,1000"
        self.area3 = "0,1000"
        self.area4 = "0,1000"
        self.area5 = "0,1000"
        self.area6 = "0,1000"
        self.area7 = "0,1000"
        self.area8 = "0,1000"
        self.area9 = "0,1000"
        self.area10 = "0,1000"
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)

    def read_colinearity(self):
        colinearity = {}
        list0 = []
        if self.colinearity_type == 'WGDI':
            num = 0
            colinearity_list = []
            for line in open(self.colinearity_file, 'r', encoding='utf-8'):
                if line[:11] == '# Alignment':
                    if len(colinearity_list) == 0:
                        pass
                    else:
                        dic_0 = {}
                        dic_0['length'] = length
                        dic_0['name'] = num
                        dic_0['pvalue'] = pvalue
                        dic_0['list'] = colinearity_list
                        colinearity[num] = dic_0
                        list0.append(num)
                        colinearity_list = []
                    num += 1
                    name = num
                    lt = line.strip('\n').split()
                    dic_ = {}
                    for i in lt:
                        if '=' in i:
                            lt0 = i.split('=')
                            dic_[str(lt0[0])] = lt0[1]
                    length = int(dic_['N'])
                    pvalue = float(dic_['pvalue'])
                else:
                    lt = line.strip('\n').split()
                    colinearity_list.append([lt[0],lt[2]])

        elif self.colinearity_type == 'MCScanX':
            num = 0
            colinearity_list = []
            for line in open(self.colinearity_file, 'r', encoding='utf-8'):
                if line[:13] == '## Alignment ':
                    # print(colinearity_list)
                    if len(colinearity_list) == 0:
                        pass
                    else:
                        dic_0 = {}
                        dic_0['length'] = length
                        dic_0['name'] = num
                        dic_0['pvalue'] = pvalue
                        dic_0['list'] = colinearity_list
                        colinearity[num] = dic_0
                        list0.append(num)
                        colinearity_list = []
                    num += 1
                    name = num
                    lt = line.strip('\n').split()
                    dic_ = {}
                    for i in lt:
                        if '=' in i:
                            lt0 = i.split('=')
                            dic_[str(lt0[0])] = lt0[1]
                    length = int(dic_['N'])
                    pvalue = float(dic_['e_value'])
                elif line[0] == '#' or line[0] == '\n':
                    continue
                else:
                    lt = line.strip('\n').split()
                    if len(lt) == 4:
                        colinearity_list.append([lt[1],lt[2]])
                    elif len(lt) == 5:
                        colinearity_list.append([lt[2],lt[3]])
                    elif len(lt) == 6:
                        colinearity_list.append([lt[3],lt[4]])
        else:
            print('blockfile_type error: File Format not recognized!')
            exit()

        if len(colinearity_list) == 0:
            pass
        else:
            dic_0 = {}
            dic_0['length'] = length
            dic_0['name'] = num
            dic_0['pvalue'] = pvalue
            dic_0['list'] = colinearity_list
            colinearity[num] = dic_0
            list0.append(num)
            colinearity_list = []


        dicx = {}
        listx = []
        for i in list0:
            dicx[i] = int(colinearity[i]['length'])
        dic1SortList = sorted(dicx.items(),key = lambda x:x[1],reverse = True)# 降序
        for i in dic1SortList:
            listx.append(i[0])
        return colinearity,listx

    def colinearity_(self,colinearity,gff,listx):# 共线性拆分
        colinearity0 = {}
        n = 0
        for i in colinearity.keys():
            lt = []
            for j in colinearity[i]['list']:
                id1,id2 = j[0],j[1]
                if gff[id1]['chr'] != gff[id2]['chr']:
                    n += 1
                    colinearity0[n] = colinearity[i]
                    break
                else:
                    if abs(gff[id1]['order']-gff[id2]['order'])<15:
                        if len(lt) == 0:
                            continue
                        else:
                            n += 1
                            colinearity0[n]=colinearity[i]
                            colinearity0[n]['length'] = len(lt)
                            colinearity0[n]['list'] = lt
                            lt = []
                    else:
                        lt.append(j)
        dicx = {}
        list0 = []
        for i in colinearity0.keys():
            dicx[i] = int(colinearity0[i]['length'])
        dic1SortList = sorted(dicx.items(),key = lambda x:x[1],reverse = True)# 降序
        for i in dic1SortList:
            # print(int(colinearity0[i[0]]['length']))
            list0.append(i[0])
        return colinearity0,list0

    def read_min_ks(self):
        if self.min_ks == '.0':# 不包括0
            op = 1
        elif self.min_ks == '0.':# 包括0
            op = 2
        else:
            op = 3


    def readks(self):
        ks = {}
        op = self.read_min_ks()
        f = open(self.ks, 'r', encoding='utf-8')
        for row in f:
            if row[0] != '#' and row[0] != '\n':
                row = row.strip('\n').split('\t')
                if row[0] != 'id1' and len(row) != 2:
                    if op == 1: 
                        if float(row[3]) > 0.0:
                            name = str(row[0]) + '_' + str(row[1])
                            ks[name] = row[3]
                        else:
                            continue
                    elif op == 2:
                        if float(row[3]) >= 0.0:
                            name = str(row[0]) + '_' + str(row[1])
                            ks[name] = row[3]
                        else:
                            continue
                    else:
                        if float(row[3]) > float(self.min_ks):
                            name = str(row[0]) + '_' + str(row[1])
                            ks[name] = row[3]
                        else:
                            continue
        return ks

    def area_(self,ks0):
        area0 = [self.area1,self.area2,self.area3,self.area4,self.area5,self.area6,self.area7,self.area8,self.area9,self.area10]
        for i in range(int(self.multiple)):
            lt = area0[i].split(',')
            if ks0 > float(lt[0]) and ks0 < float(lt[1]):
                # print(float(lt[0]),float(lt[1]))
                return True
        return False
    
    def block_ks(self,ks,lit):
        num = 0
        sum0 = 0.0
        for i in lit:
            if str(i[0]) + '_' + str(i[1]) in ks.keys():
                num += 1
                sum0 = sum0 + float(ks[str(i[0]) + '_' + str(i[1])])
            elif str(i[1]) + '_' + str(i[0]) in ks.keys():
                num += 1
                sum0 = sum0 + float(ks[str(i[1]) + '_' + str(i[0])])
            else:
                continue
        # print(sum0,num)
        if sum0 == 0:
            return 0
        avg = float(sum0/float(num))
        return avg

    def found_block(self,list0,colinearity,ks):
        out_block = []
        for i in range(len(list0)):
            if i+2 <= len(list0):
                num = 0
                sum0 = 0
                # print(colinearity[list0[i]]['list'])
                ltxx = colinearity[list0[i]]['list']+colinearity[list0[i+1]]['list']+colinearity[list0[i+2]]['list']
                for lt0 in ltxx:
                    if str(lt0[0]) + '_' + str(lt0[1]) in ks.keys():
                        ks0 = ks[str(lt0[0]) + '_' + str(lt0[1])]
                        num += 1
                        sum0 = sum0 + float(ks0)
                    elif str(lt0[1]) + '_' + str(lt0[0]) in ks.keys():
                        ks0 = ks[str(lt0[1]) + '_' + str(lt0[0])]
                        num += 1
                        sum0 = sum0 + float(ks0)
                    else:
                        pass
                if sum0 == 0:
                    continue
                ks_a = sum0/(num*1.0)
                if not self.area_(ks_a):
                    min_block_ = colinearity[list0[i+2]]['length']
                    # out_block.append(colinearity[list0[i+2]]['name'])
                    if int(min_block_) > int(self.block):
                        out_block.append(colinearity[list0[i+2]]['name'])
                        continue
                    return min_block_,out_block
            else:
                min_block_ = colinearity[list0[-1]]['length']
                return min_block_,out_block

    def writefile(self,colinearity,list0,min_block_,out_block,ks):
        # print(min_block_)
        f = open(self.savefile, 'w', encoding='utf-8')
        num = 0
        for i in range(len(list0)):
            # print(colinearity[list0[i]]['length'],list0[i])
            if int(colinearity[list0[i]]['length']) > int(min_block_) and colinearity[list0[i]]['name'] not in out_block:# and colinearity[list0[i]]['name'] not in out_block
                num += 1
                # print(num,colinearity[list0[i]]['length'],list0[i],self.block_ks(ks,colinearity[list0[i]]['list']))
                
                f.write('# Alignment '+str(num) +': pvalue='+str(colinearity[list0[i]]['pvalue'])+' N='+str(colinearity[list0[i]]['length'])+'\n')
                for n in colinearity[list0[i]]['list']:
                    f.write(str(n[0]) + '\t' + str(n[1])+'\n')
            # Alignment 1: score=259464.0 pvalue=0.0237 N=5774

    def run(self):
        gff = read_gff(self.gff)
        colinearity,list0 = self.read_colinearity()
        colinearity,list0 = self.colinearity_(colinearity,gff,list0)
        # print(colinearity)
        ks = self.readks()
        min_block_,out_block = self.found_block(list0,colinearity,ks)
        # print(out_block)
        self.writefile(colinearity,list0,min_block_,out_block,ks)
        sys.exit(0)
