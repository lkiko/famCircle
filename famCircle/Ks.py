# -*- encoding: utf-8 -*-
'''
@File        :Ks.py
@Time        :2021/09/28 11:20:37
@Author        :charles kiko
@Version        :1.0
@Contact        :charles_kiko@163.com
@Desc        :计算KS
'''


import os
import re
import random
import sys
from io import StringIO
import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO
# from Bio.Align.Applications import MafftCommandline, MuscleCommandline, ClustalwCommandline
from Bio.Phylo.PAML import yn00
from tqdm import trange
import famCircle.bez as bez
import time
import copy
# import multiprocessing # Step I : 导入模块
from multiprocessing import cpu_count#读取CPU核心数用于匹配线程数
import gc
import shutil
from multiprocessing import Pool


class Ks():
    def __init__(self, options):
        bez_conf = bez.config()
        self.pair_pep_file = 'pair.pep'
        self.pair_cds_file = 'pair.cds'
        self.prot_align_file = 'prot.aln'
        self.mrtrans = 'pair.mrtrans'
        self.pair_yn = 'pair.yn'
        self.cds_file = 'cds'
        self.pep_file = 'pep'
        self.gene_pair_linkage = 1
        self.path0 = os.getcwd()
        if cpu_count() > 32:
            self.process = 12
        else:
            self.process = cpu_count()-1
        for k, v in bez_conf:
            setattr(self, str(k), v)
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)
        self.process = int(self.process)
        self.gene_pair_linkage = int(self.gene_pair_linkage)

    def auto_file(self):
        pairs = []
        # p = pd.read_csv(self.pairs_file, sep='\t', header=None, nrows=30)
        # p = pd.read_csv(self.pairs_file, nrows=30)
        # p = '\n'.join(p[0])

        with open(self.pairs_file, 'r') as file:
            lines = [next(file) for _ in range(30)]
        p = '\n'.join(lines)

        # print(p)
        if 'path length' in p or 'MAXIMUM GAP' in p:
            collinearity = bez.read_colinearscan(self.pairs_file)
            pairs = [[v[0], v[2]] for k in collinearity for v in k[1]]
        elif 'MATCH_SIZE' in p or '## Alignment' in p:
            collinearity = bez.read_mcscanx(self.pairs_file)
            pairs = [[v[0], v[2]] for k in collinearity for v in k[1]]
        elif '# Alignment' in p:
            collinearity = bez.read_coliearity(self.pairs_file)
            pairs = [[v[0], v[2]] for k in collinearity for v in k[1]]
        elif '###' in p:
            collinearity = bez.read_jcvi(self.pairs_file)
            pairs = [[v[0], v[2]] for k in collinearity for v in k[1]]
        elif ',' in p:
            collinearity = pd.read_csv(self.pairs_file, header=None)
            pairs = collinearity.values.tolist()
        else:
            collinearity = pd.read_csv(self.pairs_file, header=None, sep='\t')
            pairs = collinearity.values.tolist()
        df = pd.DataFrame(pairs)
        df = df.drop_duplicates()
        df[0] = df[0].astype(str)
        df[1] = df[1].astype(str)
        df.index = df[0]+','+df[1]
        return df

    def run(self):
        if os.path.exists(self.ks_file):
            ks_file = open(self.ks_file, 'a+')
            ks_file.close()
        else:
            ks_file = open(self.ks_file, 'w')
            ks_file.write(
                '\t'.join(['id1', 'id2', 'ka_NG86', 'ks_NG86', 'ka_YN00', 'ks_YN00'])+'\n')
            ks_file.close()
        if os.path.exists("processes"):
            pass
        else:
            os.mkdir("processes")

        cds = SeqIO.to_dict(SeqIO.parse(self.cds_file, "fasta"))
        path = os.getcwd()
        if not os.path.exists(self.pep_file):
            if bez.cds_to_pep(os.path.join(path, self.cds_file),os.path.join(path, self.pep_file)):
                print("序列转换")
        pep = SeqIO.to_dict(SeqIO.parse(self.pep_file, "fasta"))
        print(len(pep.keys()),len(cds.keys()),"氨基酸和cds数量")
        df_pairs = self.auto_file()

        df_pairs = df_pairs[(df_pairs[0].isin(cds.keys())) & (df_pairs[1].isin(
            cds.keys())) & (df_pairs[0].isin(pep.keys())) & (df_pairs[1].isin(pep.keys()))]
        print(df_pairs)
        pairs = df_pairs[[0, 1]]#.to_numpy()# 拆分出第1、2列
        # f.write(str(len(pairs))+"\t"+"碱基对数量\n")
        # print(pairs[0][0][:3],pairs[0][1][:3])
        # if len(pairs) > 0:
        #     allpairs = []
        #     pair_hash = []
        #     print("计算价值过滤")
        #     for i in trange(len(pairs)):
        #         k = pairs[i]
        #         if k[0]+','+k[1] in pair_hash or k[1]+','+k[0] in pair_hash:
        #             continue
        #         else:
        #             pair_hash.append(k[0]+','+k[1])
        #             pair_hash.append(k[1]+','+k[0])
        #             allpairs.append(k)
        # del pairs
        print(pairs)
        print("计算价值过滤！")
        allpairs = pairs[~pairs.apply(frozenset, axis=1).duplicated()].to_numpy()#两列组合过滤
        # allpairs = random.sample(list(allpairs), len(allpairs))
        # allpairs = [list(pair) for pair in allpairs]# 
        # print(allpairs[0])
        allpairs = [tuple(pair) for pair in allpairs]# 
        # print(allpairs[0])

        allpairs = np.array(list(set(allpairs)))# 打乱
        print(allpairs)
        # exit()
        gc.collect()
        m = len(allpairs)
        print("计算量",m)
        n = int(np.ceil(m / float(self.process)))
        for i in range(int(self.process)):
            if i < int(self.process)-1:
                print(i*n,i*n+n)
            else:
                print(i*n)
        pool = Pool(int(self.process))
        for i in range(int(self.process)):
            if i < int(self.process)-1:
                new_pd = allpairs[i*n:i*n+n]
            else:
                new_pd = allpairs[i*n:]
            if os.path.exists('processes_'+str(i)):
                pass
            else:
                os.mkdir('processes_'+str(i))
            pool.apply_async(self.run0, args=(
                new_pd,i,pep,cds), error_callback=self.print_error)
        pool.close()
        pool.join()

    # 多进程错误打印
    def print_error(self,value):
        print("Process pool error, the cause of the error is :", value)

    def run0(self,pairs,n,pep,cds):
        print('The No.',n,' begin the process, process number ',os.getpid())
        os.chdir('processes_'+str(n))
        prefix = '../processes/'+str(n)
        f = open("../processes/进程详情"+str(n)+".txt",'w')
        f.write(str(len(pep.keys()))+"\t"+str(len(cds.keys()))+"\t"+"蛋白质数量 cds数量\n")
        f.write(str(len(pairs))+"\t"+"碱基对数量\n")
        # f.write(str(len(pairs))+"\t"+"碱基对数量2\n")
        # f.close()
        # for k in pairs:
        pairs = [list(tuple(pair)) for pair in pairs]
        for i in trange(len(pairs)):
            k = pairs[i]
            if k[0] == k[1]:
                continue
            for file in (prefix+self.pair_pep_file, prefix+self.pair_cds_file):
                try:
                    os.remove(file)
                except OSError:
                    pass
            SeqIO.write([cds[k[0]], cds[k[1]]], prefix+self.pair_cds_file, "fasta")
            SeqIO.write([pep[k[0]], pep[k[1]]], prefix+self.pair_pep_file, "fasta")
            if self.gene_pair_linkage > 1:# 多基因对联合
                for linkage in range(self.gene_pair_linkage-1):
                    for file in (prefix+str(linkage+1)+'_'+str(n)+'.pep', prefix+str(linkage+1)+'_'+str(n)+'.cds'):
                        try:
                            os.remove(file)
                        except OSError:
                            pass
                filtered_list = [elem for elem in pairs if elem != k]
                linkages = random.sample([i for i in range(len(filtered_list))], self.gene_pair_linkage-1)
                # print(linkages)
                for k1 in range(self.gene_pair_linkage-1):
                    k1lt = filtered_list[linkages[k1]]
                    # print(k1lt,k1)
                    SeqIO.write([cds[k1lt[0]], cds[k1lt[1]]], prefix+str(k1+1)+'_'+str(n)+".cds", "fasta")
                    SeqIO.write([pep[k1lt[0]], pep[k1lt[1]]], prefix+str(k1+1)+'_'+str(n)+".pep", "fasta")
            # exit()

            kaks = self.pair_kaks(k,prefix,n)
            # exit()
            if len(kaks) == 0:
                print(n,k,"计算错误")
                # f = open("../processes/进程详情"+str(n)+".txt",'a+')
                f.write(str(k)+"\t计算错误碱基对\n")
                # f.close()
                for file in (prefix+self.pair_pep_file, prefix+self.pair_cds_file, self.mrtrans, self.pair_yn, self.prot_align_file, '2YN.dN', '2YN.dS', '2YN.t', 'rst', 'rst1', 'yn00.ctl', 'rub'):
                    try:
                        os.remove(file)
                    except OSError:
                        pass
                continue
            ks_file = open('../'+self.ks_file, 'a+')
            ks_file.write('\t'.join([str(i) for i in list(k)+list(kaks)])+'\n')
            ks_file.close()
        # if n == int(int(self.process)) - 1:
        #     for i in trange(len(pairs)):
        #         k = pairs[i]
        #         if k[0] == k[1]:
        #             continue
        #         for file in (prefix+self.pair_pep_file, prefix+self.pair_cds_file):
        #             try:
        #                 os.remove(file)
        #             except OSError:
        #                 pass
        #         SeqIO.write([cds[k[0]], cds[k[1]]], prefix+self.pair_cds_file, "fasta")
        #         SeqIO.write([pep[k[0]], pep[k[1]]], prefix+self.pair_pep_file, "fasta")
        #         kaks = self.pair_kaks(k,prefix)
        #         if len(kaks) == 0:
        #             print(n,k,"计算错误")
        #             # f = open("../processes/进程详情"+str(n)+".txt",'a+')
        #             f.write(str(k)+"\t计算错误碱基对\n")
        #             # f.close()
        #             for file in (prefix+self.pair_pep_file, prefix+self.pair_cds_file, self.mrtrans, self.pair_yn, self.prot_align_file, '2YN.dN', '2YN.dS', '2YN.t', 'rst', 'rst1', 'yn00.ctl', 'rub'):
        #                 try:
        #                     os.remove(file)
        #                 except OSError:
        #                     pass
        #             continue
        #         ks_file = open('../'+self.ks_file, 'a+')
        #         ks_file.write('\t'.join([str(i) for i in list(k)+list(kaks)])+'\n')
        #         ks_file.close()

        # else:
        #     for k in pairs:
        #         if k[0] == k[1]:
        #             continue
        #         for file in (prefix+self.pair_pep_file, prefix+self.pair_cds_file):
        #             try:
        #                 os.remove(file)
        #             except OSError:
        #                 pass
        #         SeqIO.write([cds[k[0]], cds[k[1]]], prefix+self.pair_cds_file, "fasta")
        #         SeqIO.write([pep[k[0]], pep[k[1]]], prefix+self.pair_pep_file, "fasta")
        #         kaks = self.pair_kaks(k,prefix)
        #         if len(kaks) == 0:
        #             print(n,k,"计算错误")
        #             f.write(str(k)+"\t计算错误碱基对\n")
        #             for file in (prefix+self.pair_pep_file, prefix+self.pair_cds_file, self.mrtrans, self.pair_yn, self.prot_align_file, '2YN.dN', '2YN.dS', '2YN.t', 'rst', 'rst1', 'yn00.ctl', 'rub'):
        #                 try:
        #                     os.remove(file)
        #                 except OSError:
        #                     pass
        #             continue
        #         ks_file = open('../'+self.ks_file, 'a+')
        #         ks_file.write('\t'.join([str(i) for i in list(k)+list(kaks)])+'\n')
        #         ks_file.close()
        for file in (prefix+self.pair_pep_file, prefix+self.pair_cds_file, self.mrtrans, self.pair_yn, self.prot_align_file, '2YN.dN', '2YN.dS', '2YN.t', 'rst', 'rst1', 'yn00.ctl', 'rub'):
            try:
                os.remove(file)
            except OSError:
                pass
        gc.collect()
        f.close()
        print(n,"号测试文件写入完毕")
        
        os.chdir(self.path0)
        shutil.rmtree('processes_'+str(n))
        print('End of the No.',n,' processes, process number ',os.getpid(),'##########################################################')

    def pair_kaks(self, k,prefix,n):
        self.alignks(prefix,n)
        pal = self.pal2nal(prefix,n)
        # exit()
        if not pal:
            print("pal2nal 计算错误！")
            return []
        kaks = self.run_yn00(prefix)
        if kaks == None:
            print("yn00 计算错误！")
            return []
        kaks_new = [kaks[k[0]][k[1]]['NG86']['dN'], kaks[k[0]][k[1]]['NG86']
                    ['dS'], kaks[k[0]][k[1]]['YN00']['dN'], kaks[k[0]][k[1]]['YN00']['dS']]
        return kaks_new

    def alignks(self,prefix,n):
        if self.align_software == 'mafft':
            bez.align(self.align_software,self.mafft_path,prefix+self.pair_pep_file,self.prot_align_file)
            if self.gene_pair_linkage>1:
                for i in range(self.gene_pair_linkage-1):
                    bez.align(self.align_software,self.mafft_path,prefix+str(i+1)+'_'+str(n)+'.pep',str(i+1)+'_'+str(n)+'.pep.aln')
                prot = SeqIO.to_dict(SeqIO.parse(self.prot_align_file, "fasta"))
                seqdic = {}
                for key in prot.keys():
                    seqdic[key] = str(prot[key].seq)
                lt = sorted(list(seqdic.keys()))
                for i in range(self.gene_pair_linkage-1): 
                    prot1 = SeqIO.to_dict(SeqIO.parse(str(i+1)+'_'+str(n)+'.pep.aln', "fasta"))
                    lt1 = sorted(list(prot1.keys()))
                    # print(lt1,i,'mafft')
                    for key1,key2 in zip(lt1,lt):
                        # print()
                        seqdic[key2] = seqdic[key2]+str(prot1[key1].seq)
                f = open(self.prot_align_file,'w')
                f.write('>'+lt[0]+'\n'+seqdic[lt[0]]+'\n>'+lt[1]+'\n'+seqdic[lt[1]])
                f.close()

        if self.align_software == 'muscle':
            bez.align(self.align_software,self.muscle_path,prefix+self.pair_pep_file,self.prot_align_file)
            if self.gene_pair_linkage>1:
                for i in range(self.gene_pair_linkage-1):
                    bez.align(self.align_software,self.muscle_path,prefix+str(i+1)+'_'+str(n)+'.pep',str(i+1)+'_'+str(n)+'.pep.aln')
                prot = SeqIO.to_dict(SeqIO.parse(self.prot_align_file, "fasta"))
                seqdic = {}
                for key in prot.keys():
                    seqdic[key] = str(prot[key].seq)
                lt = sorted(list(seqdic.keys()))
                for i in range(self.gene_pair_linkage-1): 
                    prot1 = SeqIO.to_dict(SeqIO.parse(str(i+1)+'_'+str(n)+'.pep.aln', "fasta"))
                    lt1 = sorted(list(prot1.keys()))
                    for key1,key2 in zip(lt1,lt):
                        # print()
                        seqdic[key2] = seqdic[key2]+str(prot1[key1].seq)
                f = open(self.prot_align_file,'w')
                f.write('>'+lt[0]+'\n'+seqdic[lt[0]]+'\n>'+lt[1]+'\n'+seqdic[lt[1]])
                f.close()

        if self.align_software == "clustalw":
            bez.align(self.align_software,self.clustalw_path,prefix+self.pair_pep_file,self.prot_align_file)
            if self.gene_pair_linkage>1:
                for i in range(self.gene_pair_linkage-1):
                    bez.align(self.align_software,self.clustalw_path,prefix+str(i+1)+'_'+str(n)+'.pep',str(i+1)+'_'+str(n)+'.pep.aln')
                prot = SeqIO.to_dict(SeqIO.parse(self.prot_align_file, "fasta"))
                seqdic = {}
                for key in prot.keys():
                    seqdic[key] = str(prot[key].seq)
                lt = sorted(list(seqdic.keys()))
                for i in range(self.gene_pair_linkage-1): 
                    prot1 = SeqIO.to_dict(SeqIO.parse(str(i+1)+'_'+str(n)+'.pep.aln', "fasta"))
                    lt1 = sorted(list(prot1.keys()))
                    for key1,key2 in zip(lt1,lt):
                        # print()
                        seqdic[key2] = seqdic[key2]+str(prot1[key1].seq)
                f = open(self.prot_align_file,'w')
                f.write('>'+lt[0]+'\n'+seqdic[lt[0]]+'\n>'+lt[1]+'\n'+seqdic[lt[1]])
                f.close()

    def pal2nal(self,prefix,n):

        if self.gene_pair_linkage>1:
            cds = SeqIO.to_dict(SeqIO.parse(prefix+self.pair_cds_file, "fasta"))
            seqdic = {}
            for key in cds.keys():
                seqdic[key] = str(cds[key].seq)
            lt = sorted(list(seqdic.keys()))
            for i in range(self.gene_pair_linkage-1): 
                cds1 = SeqIO.to_dict(SeqIO.parse(prefix+str(i+1)+'_'+str(n)+'.cds', "fasta"))
                lt1 = sorted(list(cds1.keys()))
                # print(lt1,i,'pal2nal')
                for key1,key2 in zip(lt1,lt):
                    if seqdic[key2][-3:].upper() in ['TAA','TAG','TGA']:# 去除中间的终止密码子
                        seqdic[key2] = seqdic[key2][:-3]+str(cds1[key1].seq)
                    else:
                        seqdic[key2] = seqdic[key2][:-3]+str(cds1[key1].seq)
            f = open(prefix+self.pair_cds_file,'w')
            f.write('>'+lt[0]+'\n'+seqdic[lt[0]]+'\n>'+lt[1]+'\n'+seqdic[lt[1]])
            f.close()
        # exit()

        args = ['perl', self.pal2nal_path, self.prot_align_file,
                prefix+self.pair_cds_file, '-output paml -nogap', '>'+self.mrtrans]
        command = ' '.join(args)
        try:
            os.system(command)
            if os.path.getsize(self.mrtrans) == 0:
                print("pal2nal 计算结果为空")
                return False
        except:
            return False
        return True

    def run_yn00(self,prefix):
        yn = yn00.Yn00()
        yn.alignment = self.mrtrans
        yn.working_dir = os.getcwd()
        yn.out_file = self.pair_yn
        yn.set_options(icode=0, commonf3x4=0, weighting=0, verbose=1)
        try:
            run_result = yn.run(command=self.yn00_path)
        except:
            run_result = None
        return run_result
