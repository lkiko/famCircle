# -*- encoding: utf-8 -*-
'''
@File        :tree.py
@Time        :2022/03/09 11:17:18
@Author        :charles kiko
@Version        :1.0
@Contact        :charles_kiko@163.com
@Desc        :PhyML构树
'''

import re
import os
import sys
import gc
from math import *
import csv
from famCircle.bez import *
import pandas as pd
from Bio import Phylo
import matplotlib.pyplot as plt
plt.rc('font',family='Arial')

# -i seq_file_name
# 输入文件，phylip 格式的多序列比对结果。
# -d data_type default：nt
# 该参数的值为 nt, aa 或 generic。
# -b int
# 设置 bootstrap 次数。
# -m model
# 设置替代模型。 核酸的模型有： HKY85（默认的）, JC69, K80, F81, TN93, GTR ; 氨基酸的模型有：LG （默认的）, WAG, JTT, MtREV, Dayhoff, DCMut, RtREV, CpREV, VT, Blosum62, MtMam, HIVw, HIVb 。
# -f e,m or fA,fC,fG,fT
# 设置频率计算的方法。 e 表示使用比对结果中不同氨基酸或碱基出现的频率来计算； m 表示使用最大似然法计算碱基频率，或使用替换模型计算氨基酸频率； fA,fC,fG,fT 则是 4 个浮点数，表示 4 中碱基的频率，仅适合核酸序列。
# -v prop_invar
# 设置不变位点的比例，是一个[0,1]区间的值。或者使用 e 表示程序获得其最大似然估计值。
# -a gamma
# gamma 分布的参数。此参数值是个正数，或者使用 e 表示程序获得其最大似然估计值。在 ProtTest 软件给出的最优模型中含有 G 时，使用该参数。
# -o params
# 参数优化的选项。t 表示对 tree topology 进行优化； l 表示对 branch length 进行优化； r 表示对 rate parameters 优化。
# params=tlr 这表示对 3 者都进行优化。 params=n 表示不进行优化。

# from Bio.Phylo.Applications import PhymlCommandline
# cmd = PhymlCommandline(input='Tests/Phylip/random.phy')
# out_log, err_log = cmd()

class tree():
    def __init__(self, options):
        bez_conf = config()
        for k, v in bez_conf:
            setattr(self, str(k), v)
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)

    def phyml(self):
        phyml_exe = self.phyml_path
        command = phyml_exe + ' -i '+self.phy_file+' -d '+self.data_type+' -b '+self.bootstrap+' -m '+self.model+' -f '+self.f+' -v '+self.prop_invar+' -a '+self.gamma + ' -o '+self.params
        a = os.system(command)  #使用a接收返回值
        print(a)
        tree = Phylo.read(''.join(self.phy_file.split('.')[:-1])+".phy_phyml_tree", 'newick')
        Phylo.draw_ascii(tree)
        return tree

    def iqtree(self):
        iqtree_exe = self.iqtree_path
        command = iqtree_exe + ' -s '+self.phy_file+' -m MFP'
        a = os.system(command)  #使用a接收返回值
        print(a)
        tree = Phylo.read(''.join(self.phy_file.split('.')[:-1])+".phy.treefile", 'newick')
        Phylo.draw_ascii(tree)
        return tree

    def fasttree(self):
        fasttree_exe = self.fasttree_path
        if self.data_type == 'nt':
            command = fasttree_exe + ' -nt '+self.phy_file+' > ' + ''.join(self.phy_file.split('.')[:-1])+".phy_phyml_tree"
        else:
            command = fasttree_exe + self.phy_file+' > '+''.join(self.phy_file.split('.')[:-1])+".phy_phyml_tree"
        a = os.system(command)  #使用a接收返回值
        print(a)
        tree = Phylo.read(''.join(self.phy_file.split('.')[:-1])+".phy_phyml_tree", 'newick')
        Phylo.draw_ascii(tree)
        return tree

    def run(self):
        if self.tree_software == 'phyml':
            tree = self.phyml()
        elif self.tree_software == 'iqtree':
            tree = self.iqtree()
        elif self.tree_software == 'fasttree':
            tree = self.fasttree()
        else:
            print("famCircle 目前支持phyml/iqtree/fasttree")
            exit()
        tree.rooted=True
        Phylo.draw(tree)
        # tree = tree.as_phyloxml()
        # # tree.root.color = "gray"
        # # mcra = tree.common_ancestor({"name":"E"}, {"name":"F"})
        # # mcra.color = "salmon"
        # # tree.clade[0, 1].color = "blue"
        # # Phylo.draw(tree)
        # Phylo.write(tree, "tree.xml", "phyloxml")
        plt.savefig(self.savefile,dpi=600, bbox_inches = 'tight')