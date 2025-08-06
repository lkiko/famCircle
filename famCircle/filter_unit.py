#coding = utf-8
###
###by charles lan###
###邮箱:charles_kiko@163.com###
###

#blast数据取除自身之外的5个
import sys
from tqdm import trange
import numpy as np
import pandas as pd

class filter_unit():
    def __init__(self, options):
        self.number = 10
        self.file_type = "BLAST"
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)

    def blasts(self,df):
        dic0 = df.groupby(0).groups
        job = len(dic0)
        i = 0
        for key in dic0.keys():
            p = int(((i+1)/job)*100)
            print("\r["+"*"*p+" "*(100-p)+"] "+str(i)+"/"+str(job),end="")
            local = df.loc[dic0[key]].sort_values(by=[0,11],ascending= [True, False])
            local.reset_index(drop=True,inplace=True)
            local = local.head(int(self.number))
            local.to_csv(self.outfile, index=False, header=False,sep='\t',mode='a+')
            i += 1

        # for index,row in df.iterrows():
        #     row = list(map(str, row))
        #     if row[0] not in blast:
        #         blast.append(row[0])
        #         f.write('\t'.join(row)+'\n')
        #         number = 1
        #     else:
        #         number += 1
        #         if number <= int(self.number):
        #             f.write('\t'.join(row)+'\n')

    def read_blast0(self):# 矩阵方法
        print("读取BLAST")
        df = pd.read_csv(self.blast,header = None, sep='\t')
        df1 = df.drop(df[df[0] == df[1]].index)
        # df1.sort_values(by=[0, 11], inplace=True, ascending=[True, False])
        return df1

    def read_ks0(self):# 矩阵方法
        print("读取KS")
        df = pd.read_csv(self.blast,header = None, sep='\t')
        df1 = df.drop(df[df[0] == df[1]].index)
        df1.sort_values(by=[0, 3], inplace=True, ascending=[True, True])
        return df1

    def run(self):
        with open(self.outfile,mode="w",encoding="utf-8") as f:  #写文件,当文件不存在时,就直接创建此文件
            pass
        f.close()
        if self.file_type == 'BLAST':
            df = self.read_blast0()
        elif self.file_type == 'KS':
            df = self.read_ks0()
        print(df)
        self.blasts(df)
        # f.close()
        print()
        sys.exit(0)
