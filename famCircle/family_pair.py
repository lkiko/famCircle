# -*- encoding: utf-8 -*-
'''
@File        :family_pair.py
@Time        :2021/09/28 11:23:26
@Author        :charles kiko
@Version        :1.0
@Contact        :charles_kiko@163.com
@Desc        :基因家族关系提取
'''

import csv
import sys
from famCircle.bez import *
from tqdm import trange

class family_pair():
    def __init__(self, options):
        self.pairfile_type = 'BLAST'
        self.class1 = False
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)

    def read_family(self):
        family_list0 = []
        for i in open(self.family_list,'r',encoding='utf-8'):
            if i[0] != '\n' and i[0] != '#':
                line = i.strip('\n').split()
                if line[0] not in family_list0:
                    family_list0.append(line[0])
        return family_list0

    def run(self):
        list0 = self.read_family()
        f = open(self.savefile,'w',encoding='utf-8')
        if self.pairfile_type == 'BLAST':
            num = 0
            name_list = []
            file = open(self.gene_pair,'r',encoding='utf-8').readlines()
            for i in trange(len(file)):
                line = file[i]
                if i != '\n' and i[0] != '#':
                    line = i.strip('\n').split()
                    if line[0] == line[1]:
                        continue
                    if line[0] not in name_list:
                        name_list.append(line[0])
                        if line[0] in list0 or line[1] in list0:
                            f.write(i)
                            num = 1
                        else:
                            continue
                    else:
                        if num <= int(self.parameter):
                            if line[0] in list0 or line[1] in list0:
                                f.write(i)
                                num = 1
                            else:
                                continue
                        else:
                            continue
        if self.pairfile_type == 'MCScanX':
            file = open(self.gene_pair,'r',encoding='utf-8').readlines()
            for i in trange(len(file)):
                line = file[i]
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
                    if N >= int(self.parameter):
                        self.class1 = True
                    else :
                        self.class1 = False
                elif ('#' not in row):
                    if self.class1:
                        lt = row.strip('\n').split()
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
                        if id1 in list0 or id2 in list0:
                            f.write(id1+'\t'+id2+'\n')

        elif self.pairfile_type == 'wgdi':
            file = open(self.gene_pair,'r',encoding='utf-8').readlines()
            for i in trange(len(file)):
                line = file[i]
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
                    if N >= int(self.parameter):
                        self.class1 = True
                    else :
                        self.class1 = False
                else:
                    if self.class1:
                        lt = row.strip('\n').split()
                        if lt[0] in list0 or lt[2] in list0:
                            f.write(str(lt[0])+'\t'+str(lt[2])+'\n')
        elif self.pairfile_type == 'ColinearScan':
            file = open(self.gene_pair,'r',encoding='utf-8').readlines()
            for i in trange(len(file)):
                line = file[i]
                if (row[0] == '\n' and row[0] == '+'):
                    continue
                elif (row[:3] == 'the'):
                    lt = row.strip('\n').split()
                    N = int(lt[-1])
                    if N >= int(self.parameter):
                        self.class1 = True
                    else :
                        self.class1 = False
                else:
                    if self.class1:
                        lt = row.strip('\n').split()
                        if lt[0] in list0 or lt[2] in list0:
                            f.write(str(lt[0])+'\t'+str(lt[2])+'\n')
        else:
            print('pairfile_type error: File Format not recognized!')
            exit()
