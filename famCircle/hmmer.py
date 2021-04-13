###hmmer 用隐马尔可夫模型来对pep中的序列进行筛选保留value小于1*10-20，得到一个基因家族初等的列表，然后根据列表提取出序列，\
###再用clustalw就行序列比对,再用生成的aln文件再构建隐马尔可夫模型，用新的隐马尔可夫模型再进行筛选。再提取。
# -*- coding: UTF-8 -*-

import Bio
from glob import glob
from Bio.Align.Applications import ClustalwCommandline
from Bio import Seq
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
import subprocess
from subprocess import Popen


import re
import os
from math import *
import csv
import pandas as pd
from matplotlib.patches import *
from pylab import *
from famCircle.bez import *


class hmmer():
    def __init__(self, options):
        base_conf = config()
        for k, v in base_conf:
            setattr(self, str(k), v)
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)

    def readlist(self,filename,value):
        name = []
        f = open(filename, 'r', encoding='utf-8')
        for row in f:
            if row != '\n' and row[0] != '#':
                row = row.strip('\n').split()
                if eval(row[6]) < value:
                    if ('>' + str(row[0]) not in name):
                        name.append('>' + str(row[0]))
        return name

    def readpep(self):
        peplist = {}
        f =open(self.pep, 'r', encoding='utf-8')
        s = ''
        name = ''
        for row in f:
            if row[0] == '>':
                if s != '' and name != '':
                    peplist[name] = s
                    s = ''
                    name = ''
                name = row.strip('\n')
            elif row == '\n':
                pass
            else:
                s = s + row
        if s != '' and name != '':
            peplist[name] = s
        del name, s
        return peplist

    def writepep(self,name,seq,file):
        f = open(file, 'a+', encoding='utf-8')
        f.write(name + '\n')
        f.write(seq)
        f.close()

    def run(self):
        hmmer = '/usr/bin/'
        hmmsearch = hmmer + 'hmmsearch'
        hmmbuild = hmmer + 'hmmbuild'
        m1 = hmmsearch + ' --cut_tc --domtblout one.out ' + self.hmmmold + ' ' + self.pep
        if os.path.exists('one.out'):
            os.remove ('one.out')
        d = os.popen('sudo ' + m1).read().strip()
        # x = input()######################
        list1 = self.readlist('one.out',100)# 第一次筛选
        os.remove ('one.out')
        # print(list1)
        if len(list1) == 0:
            print('No data fit the model')
            exit()
        peplist = self.readpep()
        if os.path.exists('one.pep'):
            os.remove ('one.pep')
        for i in list1:
            if i in peplist.keys():
                seq = peplist[i]
                self.writepep(i,seq,'one.pep')
            else:
                pass
        # Clustalw = '/usr/bin/clustalw'
        in_file = 'one.pep'
        if os.path.exists('out.aln'):
            os.remove ('out.aln')
        out_file = 'out.aln'
        Clustalw_cline = ClustalwCommandline(cmd = self.clustalw_path, infile=in_file, outfile=out_file, align=True, outorder="ALIGNED", convert=True, output="pir")
        # x = input()
        a, b = Clustalw_cline()
        os.remove ('one.pep')
        os.remove ('one.dnd')
        # x = input()
        # print(Clustalw_cline)
        m2 = hmmbuild + ' ' + self.newmold + '  out.aln'
        # print(m2)
        d = os.popen('sudo ' + m2).read().strip()
        # print(m2)
        m3 = hmmsearch + ' --cut_tc --domtblout ' + self.hmmlist + ' ' + self.newmold + ' ' + self.pep
        # x = input()##################
        d = os.popen('sudo ' + m3).read().strip()
        # print(m3)
        list2 = self.readlist(self.hmmlist,1)# 第二次筛选
        print('Number of gene families: ',len(list2))
        peplist = self.readpep()
        if len(peplist) == 0:
            print('No data fit the newmodel')
            exit()
        for i in list2:
            if i in peplist.keys():
                seq = peplist[i]
                file = self.savefile
                self.writepep(i, seq, file)
            else:
                pass