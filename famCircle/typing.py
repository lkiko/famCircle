# -*- encoding: utf-8 -*-
'''
@File        :typing.py
@Time        :2021/09/28 09:03:24
@Author        :charles kiko
@Version        :1.0
@Contact        :charles_kiko@163.com
@Desc        :修改格式
'''


import os
import numpy as np
import pandas as pd

class typing():
    def __init__(self, options):
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)

    def readoutfile(self):
        out_list = str(self.domainlist).split(',')
        fam_out = []
        for i in out_list:
            pathlist = self.domainpath + '/' + str(i)
            f = open(pathlist, 'r', encoding='utf-8')
            for row in f:
                if row[0] != '#' and row[0] != '\n':
                    row = row.strip('\n').split()
                    fam_out.append(row)
        df = pd.DataFrame(fam_out)
        df.rename(columns={0: 'target_name', 1: 'accession',
                                2: 'tlen', 3: 'query_name', 4: 'accession1', 5: 'qlen',
                                6: 'E-value', 7: 'score', 8: 'bias', 9: '#', 10: 'of', 11: 'c-Evalue',
                                12: 'i-Evalue', 13: 'score2', 14: 'bias2', 15: 'from', 16: 'to',
                                17: 'from1', 18: 'to1', 19: 'from2', 20: 'to2', 21: 'acc',
                                22: 'description_of_target'}, inplace=True)
        df.drop_duplicates(subset = 'target_name')# 去重
        return df

    def run(self):
        df = self.readoutfile()
        fam_name = np.array(df['target_name'].astype(str))
        lt = []
        for i in fam_name:
            if i not in lt:
                lt.append(i)
                with open(self.savefile, "a", newline='', encoding='utf-8') as file:
                    file.write(i + '\n')
                    file.close()