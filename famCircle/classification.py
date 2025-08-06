# -*- encoding: utf-8 -*-
'''
@File        :tamdem.py
@Time        :2021/09/28 11:17:18
@Author        :charles kiko
@Version        :1.0
@Contact        :charles_kiko@163.com
@Desc        :放射型圈图
'''



import os
import sys
from math import *
from pylab import *
from famCircle.bez import *
plt.rc('font',family='Arial')


class classification():
    def __init__(self, options):
    ##### circle parameters
        self.GAP_RATIO = 5 #gaps between chromosome circle, chr:gap = 4: 1
        self.dpi = 600
        self.block= 0.008 #block scize
        self.radius = 0.3
        self.chrscize = 0.01
        self.block_scize= 0.004 #block scize
        self.blockhight = 0.005
        self.blockthick = 0.007 #0.006
        self.ks_file_type = 'ks'
        self.model = 'gradual'
        self.position = 'order'
        self.start_list = {}
        self.outfile = "classification.out"
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)

    def ksrun(self):
        lens,chr_list = read_lens(self.lens)# 染色体字典
        gff = read_gff(self.gff,chr_list)# 基因字典
        family,family_dic,dic_motif = read_family1(self.genefamily)# 家族列表
        family = drop_none(family,gff)
        colineartly = self.readblast()
        return gff,family,colineartly

    def readblast(self):
        blast_dic = {}
        one_gene = []
        alphagenepairs = open(self.genepairs, 'r', encoding='utf-8')
        if self.genepairsfile_type == 'famCircle':
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
                    if N >= int(self.block):
                        self.class1 = True
                    else :
                        self.class1 = False
                else:
                    if self.class1:
                        lt = row.strip('\n').split()
                        id1, id2 = str(lt[0]),str(lt[1])
                        if id1 not in one_gene:
                            one_gene.append(id1)
                            blast_dic[id1] = [id2]
                        else:
                            blast_dic[id1].append(id2)
        elif self.genepairsfile_type == 'WGDI':
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
                    if N >= int(self.block):
                        self.class1 = True
                    else :
                        self.class1 = False
                else:
                    if self.class1:
                        lt = row.strip('\n').split()
                        id1, id2 = str(lt[0]),str(lt[2])
                        if id1 not in one_gene:
                            one_gene.append(id1)
                            blast_dic[id1] = [id2]
                        else:
                            blast_dic[id1].append(id2)
        elif self.genepairsfile_type == 'ColinearScan':
            for row in alphagenepairs:
                if (row[0] == '\n',row[0] == '+'):
                    continue
                elif (row[:3] == 'the'):
                    lt = row.strip('\n').split()
                    N = int(lt[-1])
                    if N >= int(self.block):
                        self.class1 = True
                    else :
                        self.class1 = False
                else:
                    if self.class1:
                        lt = row.strip('\n').split()
                        id1, id2 = str(lt[0]),str(lt[2])
                        if id1 not in one_gene:
                            one_gene.append(id1)
                            blast_dic[id1] = [id2]
                        else:
                            blast_dic[id1].append(id2)
        elif self.genepairsfile_type == 'MCScanX':
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
                    if N >= int(self.block):
                        self.class1 = True
                    else :
                        self.class1 = False
                elif ('#' not in row):
                    if self.class1:
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
        elif self.genepairsfile_type == 'BLAST':
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
                        if num <= int(self.block):
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
        collinearity = []
        for key in blast_dic.keys():
            collinearity.append(key)
            for i in blast_dic[key]:
                collinearity.append(i)
        return collinearity


    def rad_to_coord(self, angle, radius):
        return radius*cos(angle), radius*sin(angle)

    def to_deg(self, bp, total):
        # from basepair return as degree
        return bp*360./total

    def to_radian(self, bp, total):
        # from basepair return as radian
        # print ("to_radian", bp, total)
        return radians(bp*360./total)

    def plot_arc(self, start, stop, radius):
        # start, stop measured in radian 染色体
        #print start, stop
        t = arange(start, stop, pi/720.)
        x, y = radius*cos(t), radius*sin(t)
        plot(x, y, "k-",alpha=.5)

    def plot_cap(self, angle, clockwise):
        radius=self.sm_radius
        # angle measured in radian, clockwise is boolean 鸭舌
        if clockwise: 
            t = arange(angle, angle+pi, pi/30.)
        else: 
            t = arange(angle, angle-pi, -pi/30.)
        x, y = radius*cos(t), radius*sin(t)
        middle_r = (self.radius_a+self.radius_b)/2
        x, y = x + middle_r*cos(angle), y + middle_r*sin(angle)
        plot(x, y, "k-", alpha=.5)

    def plot_arc_block(self, start, radius):# block
        t = arange(start, start+self.block, pi/720.)
        # print(start, start+self.block)
        x,y = radius * cos(t), radius*sin(t)
        x1, y1 = (radius-self.blockthick) * cos(t), (radius-self.blockthick) * sin(t)
        plot(x, y, "b-", linewidth=3, alpha=0.9)

    def zj(self,lens,chr_list):
        fullchrolen = sum([lens[i]['end'] for i in lens.keys()])
        fullgene = sum([lens[i]['order'] for i in lens.keys()])
        gene_average = int(fullchrolen/fullgene)
        chr_number = len(lens.keys()) # total number of chromosomes
        GAP = fullchrolen/self.GAP_RATIO/chr_number # gap size in base pair
        total_size = fullchrolen + chr_number * GAP # base pairs 
        # print('total_size',total_size)
        for i in chr_list:
            if i == chr_list[0]:
                self.start_list[i] = 0
            else:
                self.start_list[i] = self.start_list[chr_list[chr_list.index(i)-1]] + lens[chr_list[chr_list.index(i)-1]]['end'] + GAP
        return total_size,gene_average


    def transform_deg(self, ch, pos, total_size):
        return self.to_deg(pos + self.start_list[ch], total_size)

    def transform_pt(self, ch, pos, r, total_size):
        # convert chromosome position to axis coords
        rad = self.to_radian(pos + self.start_list[ch], total_size)
        return r*cos(rad), r*sin(rad)

    def transform_pt1(self,ch, pos, r, total_size):
        # convert chromosome position to axis coords
        # print("transform", ch, pos, r)
    #    print "startlist", self.start_list[ch]
        rad = self.to_radian(pos + self.start_list[ch], total_size)
        return r*cos(rad), r*sin(rad),rad

    def transform_pt2(self, rad, r):
        return r*cos(rad), r*sin(rad)

    def plot_bez_inner(self, p1, p2, cl, total_size,lw0):
    #    print "inner"
        a, b, c = p1
        ex1x, ex1y = self.transform_pt(a, b, c, total_size)
        a, b, c = p2
        ex2x, ex2y = self.transform_pt(a, b, c, total_size)
        # Bezier ratio, controls curve, lower ratio => closer to center
        ratio = .5
        x = [ex1x, ex1x*ratio, ex2x*ratio, ex2x]
        y = [ex1y, ex1y*ratio, ex2y*ratio, ex2y]
        step = .01
        t = arange(0, 1+step, step)
        xt = Bezier(x, t)# 贝塞尔曲线
        yt = Bezier(y, t)
        plot(xt, yt, '-', color=cl, lw=lw0, alpha=0.3)#alpha 透明度

    def plot_bez_Ks2(self, rad1, r1, rad2, r2, col, ratio):
    #    print "bez Ks 2"
        # print('rad1, r1',rad1, r1)
        ex1x, ex1y = self.transform_pt2(rad1, r1)
        ex2x, ex2y = self.transform_pt2(rad2, r2)
        # ratio = -0.7#0.5
        sita = pi / 2
        if ex1x != ex2x:
            sita = atan((ex2y-ex1y)/(ex2x-ex1x))
        d = sqrt((ex2x-ex1x)**2+(ex2y-ex1y)**2)
        L = d * ratio
        P1x = ex1x + L*sin(sita)
        P1y = ex1y - L*cos(sita)
        P2x = ex2x + L*sin(sita)
        P2y = ex2y - L*cos(sita)
        step = .01
        t = arange(0, 1+step, step)
        x=[ex1x, P1x, P2x, ex2x]
        y=[ex1y, P1y, P2y, ex2y]
        # print('x,y,t',x,y,t)
        xt = Bezier(x,t)
        yt = Bezier(y,t)
        plot(xt, yt, '-', color = col, lw = 0.7)#0.1

    def cluster0(self,cluster,id1,id2):
        for i in cluster:
            if id1 in i and id2 in i:
                return True
            else:
                continue
        return False

    def run(self):
        if (os.path.exists(self.outfile)):
            os.remove(self.outfile)
        gff,family,colineartly = self.ksrun()
        print("家族成员分类")
        print(len(family))
        WGD_segmental = []
        tandem = []
        proximal = []
        dispersed = []
        classification_dic = {}


        ### tandem and WGD监测
        fam_dic = {}
        for i in family:
            if i in colineartly:
                WGD_segmental.append(i)
                classification_dic[i] = 'WGD_segmental'
                continue
            if gff[i]['chr'] not in fam_dic.keys():
                dic0 = {}
                dic0[gff[i]['order']] = i
                fam_dic[gff[i]['chr']] = dic0
            else:
                fam_dic[gff[i]['chr']][gff[i]['order']] = i

        for i in fam_dic.keys():
            lt = list(fam_dic[i].keys())
            lt.sort()# 升序
            cluster0 = []
            # print(lt)
            for j in lt:
                if j+1 in lt:
                    tandem.append(fam_dic[i][j])
                    tandem.append(fam_dic[i][j+1])
                    classification_dic[fam_dic[i][j]] = 'tandem'
                    classification_dic[fam_dic[i][j+1]] = 'tandem'
                elif j-1 in lt:
                    tandem.append(fam_dic[i][j])
                    tandem.append(fam_dic[i][j-1])
                    classification_dic[fam_dic[i][j]] = 'tandem'
                    classification_dic[fam_dic[i][j-1]] = 'tandem'
                else:
                    cluster0.append(fam_dic[i][j])
            fam_dic[i] = cluster0
        tandem = list(set(tandem))
        print(len(tandem),len(WGD_segmental))
        # print(fam_dic)
        print("111111111")


        ####proximal检测
        fam_dic1 = {}
        for chr0 in fam_dic.keys():
            for gene in fam_dic[chr0]:
                if chr0 not in fam_dic1.keys():
                    dic0 = {}
                    dic0[gff[gene]['order']] = gene
                    fam_dic1[chr0] = dic0
                else:
                    fam_dic1[chr0][gff[gene]['order']] = gene

        # print(fam_dic1)
        fam_dic = fam_dic1
        print('2222222222')

        for i in fam_dic.keys():
            lt = list(fam_dic[i].keys())
            lt.sort()# 升序
            cluster0 = []
            # print(lt)
            for j in lt:
                for n in range(2,20):
                    if j+n in lt:
                        proximal.append(fam_dic[i][j])
                        proximal.append(fam_dic[i][j+n])
                        classification_dic[fam_dic[i][j]] = 'proximal'
                        classification_dic[fam_dic[i][j+n]] = 'proximal'
                    elif j-n in lt:
                        proximal.append(fam_dic[i][j])
                        proximal.append(fam_dic[i][j-n])
                        classification_dic[fam_dic[i][j]] = 'proximal'
                        classification_dic[fam_dic[i][j-n]] = 'proximal'
        proximal = list(set(proximal))
        print(len(proximal))
        print("3333333")

        for i in family:
            if i not in tandem and i not in proximal and i not in WGD_segmental:
                dispersed.append(i)
                classification_dic[i] = 'dispersed'
        print(len(dispersed))
        print(len(classification_dic),classification_dic)
        out = open(self.outfile,'w')
        for gene in family:
            lt = [gene,classification_dic[gene]]
            out.write('\t'.join(lt)+'\n')
        out.close()
        # sys.exit(0)
