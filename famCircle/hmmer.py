###hmmer 用隐马尔可夫模型来对pep中的序列进行筛选保留value小于1*10-20，得到一个基因家族初等的列表，然后根据列表提取出序列，\
###再用clustalw就行序列比对,再用生成的aln文件再构建隐马尔可夫模型，用新的隐马尔可夫模型再进行筛选。再提取。
# python
# -*- encoding: utf-8 -*-
'''
@File        :hmmer.py
@Time        :2021/09/28 11:23:06
@Author        :charles kiko
@Version        :1.0
@Contact        :charles_kiko@163.com
@Desc        :扫描hmm
'''



import Bio
from glob import glob
from io import StringIO
# from Bio.Align.Applications import MafftCommandline, MuscleCommandline, ClustalwCommandline
from Bio import Seq
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
import subprocess
from subprocess import Popen
import shutil
import glob # 都是标准库的东西

import re
import os
from math import *
import csv
import pandas as pd
from matplotlib.patches import *
from pylab import *
from famCircle.bez import *
#system 阻塞
import platform
import time

class hmmer():
    def __init__(self, options):
        self.cds = 'None'
        base_conf = config()
        for k, v in base_conf:
            setattr(self, str(k), v)
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)
        # self.hmmer = self.hmmer_path
        self.hmmsearch = self.hmmer_path + 'hmmsearch'
        self.hmmbuild = self.hmmer_path + 'hmmbuild'
        self.maker = '/'
        if platform.system() == 'Windows':
            self.hmmsearch = self.hmmsearch.replace("/", "\\")
            self.hmmbuild = self.hmmbuild.replace("/", "\\")
            self.maker = "\\"


    def readlist(self,filename,value):
        name = []
        f = open(filename, 'r', encoding='utf-8')
        for row in f:
            if row != '\n' and row[0] != '#':
                row = row.strip('\n').split()
                if eval(row[6]) < float(value):
                    if (str(row[0]) not in name):
                        name.append(str(row[0]))
        if len(name) == 0:
            print("Warning: Screening is too strict!")
        return name

    def readseq(self,file):
        pep = SeqIO.to_dict(SeqIO.parse(file, "fasta"))# 提取之后直接返回字典
        return pep

    def changefile(self,i):
        name = str(i) + '.new.hmm'
        oldname = str(i)
        d = os.system(self.hmmbuild+" %s %s > fromat-buile.log" % (name, oldname))
        if os.path.exists('fromat-buile.log'):
            os.remove ('fromat-buile.log')
        return name

    def runhmm(self, hmmmold):
        name = hmmmold.split(self.maker)[-1]
        os.makedirs('temp', exist_ok=True)

        if self.format_conversion == 'True':
            try:
                hmmmold = self.changefile(hmmmold)
            except:
                pass

        m1 = self.hmmsearch + ' --domtblout ./temp/one.out -E ' + self.e_value1 + ' ' + hmmmold + ' ' + self.pep + '> hmmsearch1.log'
        d = os.system(m1)# 第一次搜索
        if os.path.exists('hmmsearch1.log'):
            os.remove ('hmmsearch1.log')
        # exit()
        if os.path.exists('./temp/one.out'):
            pass
        else:
            return 0
        list1 = self.readlist('./temp/one.out',self.e_value1)# 第一次筛选
        os.remove ('./temp/one.out')
        if len(list1) == 0:
            print('STEP 1 :No data fit the model')
            return 0
        # 第一次搜索结束
        peplist = self.readseq(self.pep)
        if os.path.exists('./temp/one.pep'):
            os.remove ('./temp/one.pep')
        list1=list(set(list1))
        for i in list1:# 第一次生成筛选之后的pep
            if i in peplist.keys():
                seq = peplist[i].seq
                write_fasta(seq,i,'./temp/one.pep')
            else:
                print('蛋白质 ',i ,' 不在输入的pep文件中。')

        if self.align_software == 'mafft':
            align(self.align_software,self.mafft_path,'./temp/one.pep','./temp/one.aln')
        elif self.align_software == 'muscle':
            align(self.align_software,self.muscle_path,'./temp/one.pep','./temp/one.aln')
        elif self.align_software == "clustalw":
            align(self.align_software,self.clustalw_path,'./temp/one.pep','./temp/one.aln')

        m2 = self.hmmbuild + ' ./out_hmm/' + name + '.my.hmm ./temp/one.aln > fromat-buile2.log'
        d = os.system(m2)# 格式转换
        if os.path.exists('fromat-buile2.log'):
            os.remove ('fromat-buile2.log')
        m3 = self.hmmsearch + ' --domtblout ./out_list/' + name + '.out -E ' + self.e_value2 + ' ./out_hmm/' + name + '.my.hmm ' + self.pep + ' > hmmsearch2.log'
        d = os.system(m3)# 第二次搜索，生成out文件
        if os.path.exists('hmmsearch2.log'):
            os.remove ('hmmsearch2.log')
        hmmlist = './out_list/' + name + '.out'
        list2 = self.readlist(hmmlist,self.e_value2)# 第二次筛选
        family = open('./out_list/' + name + '.list','w')
        print('Number of gene families: ',len(list2))
        peplist = self.readseq(self.pep)
        list22 = list(set(list2))
        if len(list22) == 0:
            print('SETP 2 :No data fit the newmodel')
            return 0
        if self.cds != 'None':
            cdslist = self.readseq(self.cds)
        for i in list22:
            family.write(i+'\n')
            seq = peplist[i].seq
            file = './out_pep/'+name + '.pep'
            write_fasta(seq,i,file)
            try:
                cdsseq = cdslist[i].seq
                file = './out_cds/'+name + '.cds'
                write_fasta(cdsseq,i,file)
            except:
                print("Warning: no cdsfile or seq not in cdsfile")
        family.close()
        return 1

    def run(self):
        for file in ["out_hmm","out_cds","out_pep","out_list"]:
            os.makedirs(file, exist_ok=True)
        hmmlistname = self.hmmmoldpath.split(',')

        if self.format_conversion == 'True':
            hmmlistname1 = []
            for file in hmmlistname:
                if self.align_software == 'mafft':
                    align(self.align_software,self.mafft_path,file,file+'.aln')
                elif self.align_software == 'muscle':
                    align(self.align_software,self.muscle_path,file,file+'.aln')
                elif self.align_software == "clustalw":
                    align(self.align_software,self.clustalw_path,file,file+'.aln')
                hmmlistname1.append(file+'.aln')
            hmmlistname = hmmlistname1
        print(hmmlistname)
        for i in hmmlistname:# 遍历模型文件
            if not self.runhmm(i):
                print(i,' not find families!')
                continue
            else:
                print(i,' end!')
        folder_path = 'temp'
        try:
            shutil.rmtree(folder_path)
            print(f"{folder_path} 文件夹已成功删除。")
        except FileNotFoundError:
            print(f"{folder_path} 文件夹未找到，无法删除。")
        except Exception as e:
            print(f"删除 {folder_path} 文件夹时出现错误： {e}")
