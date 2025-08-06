# -*- encoding: utf-8 -*-
'''
@File        :msa.py
@Time        :2022/03/09 11:17:18
@Author        :charles kiko
@Version        :1.0
@Contact        :charles_kiko@163.com
@Desc        :多序列比对
'''

import re
import os
import sys
import gc
from math import *
import csv
from famCircle.bez import *
import pandas as pd

# from Bio.Align.Applications import ClustalwCommandline
# from Bio.Align.Applications import MuscleCommandline

from Bio import AlignIO

from Bio import Phylo
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class msa():
    def __init__(self, options):
    ##### circle parameters
        bez_conf = config()
        self.seq_type = 'dna'
        for k, v in bez_conf:
            setattr(self, str(k), v)
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)

# import subprocess
# from Bio import AlignIO

#     # def align(self, software, exe_path, input_file, output_file):
#     #     if software == "muscle":
#     #         cmd = [
#     #             exe_path,
#     #             "-in", input_file,
#     #             "-out", output_file,
#     #             "-clwstrict",  # 输出 clustal 格式
#     #             "-seqtype", self.seq_type  # 你自己传的是 "protein" 或 "dna"
#     #         ]
#     #     elif software == "clustalw":
#     #         cmd = [
#     #             exe_path,
#     #             "-infile=" + input_file,
#     #             "-outfile=" + output_file,
#     #             "-type=" + self.seq_type
#     #         ]
#     #     else:
#     #         raise ValueError(f"Unsupported alignment software: {software}")

#     #     # 执行命令
#     #     result = subprocess.run(cmd, capture_output=True, text=True)

#     #     if result.returncode != 0:
#     #         print("Error running", software)
#     #         print(result.stderr)
#     #         raise RuntimeError(result.stderr)

#     #     # 读取输出的clustal格式文件并保存为phylip格式（可选）
#     #     try:
#     #         alignment = AlignIO.read(output_file, "clustal")
#     #         AlignIO.write([alignment], output_file, "phylip")
#     #     except Exception as e:
#     #         print(f"Failed to parse alignment output: {output_file}")
#     #         print(e)


#     def clustalw(self):
#         clustalw_exe = self.clustalw_path
#         clustalw_cline = ClustalwCommandline(clustalw_exe, infile=self.fasta, outfile=self.outfile, type= self.seq_type)
#         stdout, stderr = clustalw_cline()
#         alignment = AlignIO.read(self.outfile, "clustal")
#         AlignIO.write([alignment], self.outfile, "phylip")
#         # align = AlignIO.read(self.outfile, "clustal")
#         # print(align)

#         # tree = Phylo.read(''.join(self.fasta.split('.')[:-1])+".dnd", "newick")
#         # Phylo.draw_ascii(tree)

#     def muscle(self):
#         muscle_exe = self.muscle_path
#         muscle_cline = MuscleCommandline(
#             cmd=muscle_exe, input=self.fasta, out=self.outfile, seqtype=self.seq_type, clwstrict=True)# dna protein
#         # print(muscle_cline)
#         stdout, stderr = muscle_cline()
#         alignment = AlignIO.read(self.outfile, "clustal")
#         AlignIO.write([alignment], self.outfile, "phylip")
#         # print(muscle_cline)
#         # align = AlignIO.read(self.outfile, "muscle")
#         # print(align)

    def run(self):
        orchid_dict = SeqIO.to_dict(SeqIO.parse(self.sequences_file, "fasta"))# 提取之后直接返回字典
        # print(orchid_dict.keys())
        data = pd.read_csv(self.clusters_file,encoding='utf-8',header=0,index_col=0,error_bad_lines=False)
        clusters = str(data.loc[int(self.number), 'gene_list']).split(',')
        my_records = []
        for i in clusters:
            rec1 = orchid_dict[i]
            my_records.append(rec1)
        SeqIO.write(my_records, self.fasta, "fasta")

        # muscle 方式
        if self.align_software == 'muscle':
            align(self.align_software, '', self.fasta, self.outfile)
            # 下面这段代码你已写得不错，可保留

        # clustalw 方式
        if self.align_software == "clustalw":
            align(self.align_software, '', self.fasta, self.outfile)
            # 同样保留你已有的多文件合并逻辑


        # if self.align_software == "clustalw":
        #     self.clustalw()
        # elif self.align_software == "muscle":
        #     self.muscle()
        # else:
        #     print("famCircle 目前支持clustalw和muscle")
