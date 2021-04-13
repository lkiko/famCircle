### draw chromosome homology famCircle for eudicots with VV as reference, if VV is provided. If no vv chro is involved, it draws famCircle for other species
# -*- coding: UTF-8 -*-
import re
import os
from math import *
import csv
import pandas as pd
from matplotlib.patches import *
from pylab import *

from famCircle.bez import *


class inner():
    def __init__(self, options):
        self.genessorted = []
        self.geneno = 0
        self.genepair2Ks = {}
        self.genepair2Ka = {}
        
        ##### circle parameters
        self.GAP_RATIO = 4 #gaps between chromosome circle, chr:gap = 4: 1
        self.radius_a, self.radius_b = .33, .335
        self.sm_radius=(self.radius_b-self.radius_a)/2 #telomere capping
        # figurefile = "_".join(sys.argv[1:len(sys.argv)])+".genefam"

        #### gene block parameters
        self.blocklength = 0.01 # radian 0.01
        self.blockthick = 0.007 #0.006
        self.shiftratio = 1.1 # define the distance between overlapping glocks
        # chrolist = sys.argv[1:len(sys.argv)]
        self.specieslist = []
        self.iscompletegenome = {}
        self.chro2len = {}
        self.vvchrolist = []
        self.otherchrolist = []
        self.labels = []
        self.start_list = []
        self.gene2pos={}
        self.gene2chain = {}
        self.peripheral = False
        self.genes = []
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)
    def readname(self, name):
        if ('^' in name):
            return name
        elif ('g' in name):
            name = name.replace('g', '^')
            return name
    # def ks(self, figurefile, fpchrolen, fpgff, fpgenefamilyinf, chrolist):
    def ksrun(self):
        # print('kaishi')
        figure(1, (8, 8))  ### define the a square, or other rectangle of the figure, if to produce an oval here
        root =axes([0, 0, 1, 1])
        fpchrolen = open(self.lens, 'r', encoding='utf-8')
        fpgff = open(self.gff, 'r', encoding='utf-8')
        fpgenefamilyinf = open(self.ks, 'r', encoding='utf-8')
        genefamily = open(self.genefamily, 'r', encoding='utf-8')
        #print "Usage: python draw.chro.famCircle.py chrolist\n"

        vvchro2color = {"VV1": "c", "VV2": "m", "VV3": "lightsalmon", "VV4": "aquamarine",
                        "VV5": "orange", "VV6": "#FFFF33", "VV7": "#FFCC33","VV8": "#FF6666",
                        "VV9": "#FF3333","VV10": "#FF3399","VV11": "#FF99FF","VV12": "#CCFFCC",
                        "VV13": "#CCFF33","VV14": "#CCCCCC","VV15": "#CCCC99","VV16": "#CC99FF",
                        "VV17": "#CC9966","VV18": "#CC9900","VV19": "#CC66FF"}

        #### self.specieslist and initial chrolist
        # print("lens")
        chrolist = []
        for row in fpchrolen:
            chro = row.split('\t')[0]
            # print(chro)
            chrolist.append(chro)
        fpchrolen.close()
        # chrolist = self.chrolist.split(',')
        for i in range(len(chrolist)):
            string = chrolist[i]
            # print (string[0:2])
            isnew = 1
            for sp in self.specieslist:
                if sp == string[0:2]:
                    isnew = 0
                    break
            if isnew==1:
                self.specieslist.append(string[0:2])
                if string == string[0:2]:
                    self.iscompletegenome[string[0:2]] = 1
                else:
                    self.iscompletegenome[string[0:2]] = 0

        ### input chromosome length
        # print("lens")
        fpchrolen = open(self.lens, 'r', encoding='utf-8')
        # print(self.iscompletegenome.keys())
        for row in fpchrolen:
            # print('ceshi')
            chro,length = row.split('\t')[0],row.split('\t')[1]
            # print(chro,'***********************')
            # chro = chro.upper()
            if len(chro) > 10 :####限制chr长度
                continue
            # sp = chro[:2].upper()
            sp = chro[:2]
            if self.iscompletegenome[sp] == 1 :
                # print(length)
                self.chro2len[chro] = int(length)
                if(re.search("VV", chro)):
                    self.vvchrolist.append(chro)
                else:
                    self.otherchrolist.append(chro)
            else:
                if chro in chrolist :
                    self.chro2len[chro] = int(length)
                    if(re.search("VV", chro)):
                        self.vvchrolist.append(chro)
                    else:
                        self.otherchrolist.append(chro)

        fpchrolen.close()
        ### full chro list
        for i in self.vvchrolist:
            self.labels.append(i)

        for i in self.otherchrolist:
            self.labels.append(i)

        # self.labels = self.vvchrolist + self.otherchrolist
        # print(self.labels)
        #print "all chro", self.labels

        ## input gff
        # gene2pos={}
        # gene2chain = {}
        familylist = []
        for row in genefamily:
            if (row != '\n'):
                familyname = row.strip('\n')
                familylist.append(self.readname(familyname))

        for row in fpgff:
            ch, gene, start, end, stran = row.split()[0],row.split()[1],row.split()[2],row.split()[3],row.split()[4]    #mrna, , software
            # ch = gene.split("^")[0]
            gene = self.readname(gene)
            if (gene in familylist):
                # print(gene)
                start = int(start)
                end = int(end)
                self.gene2pos[gene] = int(start)
                self.gene2chain[gene] = int((end-start)/abs(end -start))
                #print "+++++++++++++", gene, int(start), self.gene2chain,ch
                # print("my chro number is", ch, start)
            else:
                pass
        fpgff.close()
        
        ### input gene family gene1 gene2 Ka Ks
        # genes = []
        for row in fpgenefamilyinf:
            if (row.split()[0] == 'id1'):
                continue
            gene1, gene2, Ka, Ks = row.split()[0],row.split()[1],row.split()[2],row.split()[3]
            gene1 = self.readname(gene1)
            gene2 = self.readname(gene2)
            if (gene1 not in self.gene2pos.keys() or gene2 not in self.gene2pos.keys()):
                continue
            # print(gene1)
            genepair = ""
            # print(type(gene1))
            if gene1 < gene2:
               genepair = gene1 + " " + gene2
            else:
               genepair = gene2 + " " + gene1
            if len(gene1.split("^")[0]) < 10 and len(gene2.split("^")[0]) < 10:#数字筛选chr名字长度
               self.genepair2Ka[genepair] = float(Ka)
               self.genepair2Ks[genepair] = float(Ks)
            # print("+++++++++++++++++++++++++++++++++")
            if gene1 not in self.genes:
                if len(gene1.split("^")[0]) < 10:
                    self.genes.append(gene1)
                    # print(gene1)
            if gene2 not in self.genes :
                if len(gene2.split("^")[0]) < 10:
                    self.genes.append(gene2)
                    # print(gene1)
        fpgenefamilyinf.close()

        # print("genes ", self.genes, "genes")
        return root

    def rad_to_coord(self,angle, radius):
        return radius*cos(angle), radius*sin(angle)

    def to_deg(self,bp, total):
        # from basepair return as degree 从碱基对返回为学位
        return bp*360./total

    def to_radian(self,bp, total):# 返回弧度
        # from basepair return as radian
        # print ("to_radian", bp, total)
        return radians(bp*360./total)

    def plot_arc(self,start, stop, radius):
        # start, stop measured in radian
        #print start, stop
        t = arange(start, stop, pi/720.)
        x, y = radius*cos(t), radius*sin(t)
        plot(x, y, "k-", alpha=.5)

    def plot_cap(self,angle, clockwise=True):
        radius=self.sm_radius
        # angle measured in radian, clockwise is boolean
        if clockwise: t = arange(angle, angle+pi, pi/30.)
        else: t = arange(angle, angle-pi, -pi/30.)
        x, y = radius*cos(t), radius*sin(t)
        middle_r = (self.radius_a+self.radius_b)/2
        x, y = x + middle_r*cos(angle), y + middle_r*sin(angle)
        plot(x, y, "k-", alpha=.5)

    def plot_arc_block(self,start, chain, radius):
        t = arange(start, start+self.blocklength, pi/720.)
        x,y = radius * cos(t), radius*sin(t)
        x1, y1 = (radius-self.blockthick) * cos(t), (radius-self.blockthick) * sin(t)
        plot(x, y, "b-", alpha=0.5)
    #    plot(x1, y1, "g-", alpha=0.5)

        x0, y0 = radius*cos(start), radius*sin(start)
        x1,y1  = (radius-self.blockthick)*cos(start), (radius-self.blockthick)*sin(start)
    #    plot([x0, y0], [x1, y1], "g-", lw=0.2)

        x0, y0 = radius*cos(start+self.blocklength), radius*sin(start+self.blocklength)
        x1,y1  = (radius-self.blockthick)*cos(start+self.blocklength), (radius-self.blockthick)*sin(start+self.blocklength)
    #    plot([x0, y0], [x1, y1], "g-", lw=0.2)
    def zj(self):
        # print(self.chro2len.values())
        # print(pd.DataFrame(self.chro2len.values()).sum())
        fullchrolen = int(pd.DataFrame(self.chro2len.values()).sum())
        chr_number = len(self.labels) # total number of chromosomes
        GAP = fullchrolen/self.GAP_RATIO/chr_number # gap size in base pair
        total_size = fullchrolen + chr_number * GAP # base pairs
        
        for i in range(chr_number):
            self.start_list.append(0)
        # self.start_list = [0]*chr_number
        for i in range(1, chr_number):
            self.start_list[i] = self.start_list[i-1] + self.chro2len[self.labels[i-1]] + GAP
        stop_list = [(self.start_list[i] + self.chro2len[self.labels[i]]) for i in range(chr_number)]
        # print("start list", self.start_list)
        return stop_list, total_size, chr_number

    def transform_deg(self,ch, pos):

        return self.to_deg(pos + self.start_list[ch], total_size)

    def transform_pt(self,ch, pos, r, total_size):
        # convert chromosome position to axis coords
        # print("transform", ch, pos, r)
    #    print "startlist", self.start_list[ch]
        rad = self.to_radian(pos + self.start_list[ch], total_size)
        return r*cos(rad), r*sin(rad)
    def transform_pt2(self,rad, r):
        return r*cos(rad), r*sin(rad)

    def plot_bez_inner(self,p1, p2, cl, total_size):
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
        xt = Bezier(x, t)
        yt = Bezier(y, t)
        plot(xt, yt, '-', color=cl, lw=.1)

    def plot_bez_outer(self,p1, p2, cl):
    #    print "outer"
        a, b, c = p1
        ex1x, ex1y = self.transform_pt(a, b, c, total_size)
        a, b, c = p2
        ex2x, ex2y = self.transform_pt(a, b, c, total_size)
        # Bezier ratio, controls curve, lower ratio => closer to center
        ratio = 1.1
        x = [ex1x, ex1x*ratio, ex2x*ratio, ex2x]
        y = [ex1y, ex1y*ratio, ex2y*ratio, ex2y]
        step = .01
        t = arange(0, 1+step, step)
        xt = Bezier(x, t)
        yt = Bezier(y, t)
        plot(xt, yt, '-', color=cl, lw=.1)

    def plot_bez_Ks(self,rad1, r1, rad2, r2, col, makers):
    #    print "bez Ks 1"
        ex1x, ex1y = self.transform_pt2(rad1, r1)
        ex2x, ex2y = self.transform_pt2(rad2, r2)
        ratio = 0.5
        x = [ex1x, ex1x*ratio, ex2x*ratio, ex2x]
        y = [ex1y, ex1y*ratio, ex2y*ratio, ex2y]
        step = .01
        t = arange(0, 1+step, step)
        xt = Bezier(x, t)
        yt = Bezier(y, t)
        plot(xt, yt, makers, color=col, lw=.1)

    def plot_bez_Ks2(self,rad1, r1, rad2, r2, col, makers):
    #    print "bez Ks 2"
        ex1x, ex1y = self.transform_pt2(rad1, r1)
        ex2x, ex2y = self.transform_pt2(rad2, r2)
        ratio = 0.5
        sita = pi/2
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
        xt = Bezier(x,t)
        yt = Bezier(y,t)
        plot(xt, yt, makers, color = col, lw = 0.1)
    def plot_bez_Ks3(self,rad1, r1, rad2, r2, col, makers):
    #    print "bez Ks 3"
        ex1x, ex1y = self.transform_pt2(rad1, r1)
        ex2x, ex2y = self.transform_pt2(rad2, r2)
        ratio = -0.5
        sita = pi/2
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
        xt = Bezier(x,t)
        yt = Bezier(y,t)
        plot(xt, yt, makers, color = col, lw = 0.1)

    def plot_sector(self, s, root):
        # block_id, chr_id, start_bp, stop_bp
        block_id, chr_id, start_bp, stop_bp = s
        # self.zj()
        theta0, theta1 = self.transform_deg(chr_id, start_bp), self.transform_deg(chr_id, stop_bp)
        dtheta = .1
        p1 = Wedge((0, 0), self.radius_a, theta0, theta1, dtheta=dtheta, fc="w", ec="w")
        p2 = Wedge((0, 0), self.radius_b, theta0, theta1, dtheta=dtheta, fc=colors[block_id], ec="w")
        root.add_patch(p2)
        root.add_patch(p1)
        return root

    def run_plot(self, Ks, start1, start2):
        colorlist = ['red', 'chocolate', 'orange', 'yellow', 'palegreen', 'lime', 'teal', 'dodgerblue', 'blueviolet', 'fuchsia']
        # makerlist = ['-', '.', '-.']
        lt = self.ks_concern.strip('\n').split(',')
        length = len(lt)
        if (Ks > eval(lt[-1]) or Ks < eval(lt[0])):
            color = 'gray'
            class1 = 'None'
            if (eval(self.peripheral) == True):
                if (abs(start1 - start2) < eval(self.bridge)):
                    makers = '-'
                    return color, makers, class1
                else:
                    return None,None,None
            else:
                return None,None,None
        else:
            for i in range(length - 1):
                # min = eval(lt[i])
                max = eval(lt[i + 1])
                if (Ks <= max):
                    color = colorlist[i]
                    class1 = i
                    if (abs(start1 - start2) < eval(self.bridge)):
                        makers = '-'
                        return color, makers, class1
                    else:
                        return None,None,None

                else:
                    continue



    def run(self):
        if (os.path.exists(self.savecsv)):
            os.remove(self.savecsv)
        lt = ['gene1', 'gene2', 'start1', 'shift1', 'start2', 'shift2', 'class1']
        with open(self.savecsv, "a", newline='', encoding='utf-8') as file:
            writer = csv.writer(file ,delimiter=',')
            writer.writerow(lt)
        root = self.ksrun()
        stop_list, total_size, chr_number = self.zj()
        gene2pos = self.gene2pos
        ## sort gene according to lacation on circle
        #
        # print(len(self.genes),'++++++++++++++++++++++++')
        for i in range(len(self.genes)):
            # print(i,self.genes[i])
            # print(i,'**************************')
            # print(self.genes[i])
            if self.geneno == 0:
                
                self.genessorted.append(self.genes[0])
                self.geneno = self.geneno + 1
            else:
                firstgene = self.genessorted[0]
                # print(firstgene)
                lastgene = self.genessorted[-1]
                chf = firstgene.split("^")[0]
                chl = lastgene.split("^")[0]
             
                # print(firstgene, lastgene, chf, chl, gene2pos[firstgene])
                posf = gene2pos[firstgene] + self.start_list[self.labels.index(chf)]
                posl = gene2pos[lastgene] + self.start_list[self.labels.index(chl)]
                chi = self.genes[i].split("^")[0]
                posi = gene2pos[self.genes[i]] + self.start_list[self.labels.index(chi)]
        #       print posf, posl, posi
                if posi <= posf:
                   self.genessorted[0:0] = [self.genes[i]]
                elif posi >= posl:
                    self.genessorted.append(self.genes[i])
                else:
                    for j in range(len(self.genessorted)-1):
                        chj = self.genessorted[j].split("^")[0]
                        posj = gene2pos[self.genessorted[j]] + self.start_list[self.labels.index(chj)]
                        chj1 = self.genessorted[j+1].split("^")[0]
                        posj1 = gene2pos[self.genessorted[j+1]]+self.start_list[self.labels.index(chj1)]
                        # print(posj, posj1, posi)
                        if posi > posj and posi < posj1:
                            self.genessorted[j+1:j+1] = [self.genes[i]]
                # print("genesort ", self.genessorted)

        #print "genesorted:", self.genessorted, "self.genessorted"

        # the chromosome layout
        j = 0
        for start, stop in zip(self.start_list, stop_list):
            start, stop = self.to_radian(start, total_size), self.to_radian(stop, total_size)
            # shaft
            self.plot_arc(start, stop, self.radius_a)
            self.plot_arc(start, stop, self.radius_b)

            # telemere capping
            self.plot_cap(start, clockwise=False)
            self.plot_cap(stop)
            # chromosome self.labels
            label_x, label_y = self.rad_to_coord((start+stop)/2, self.radius_b*1.2)
            text(label_x, label_y, self.labels[j], horizontalalignment="center", verticalalignment="center")
            j+=1


        ### define shift level and draw gene blocks
        shiftlevel = 0
        gene2shift = {}
        gene2location = {}
        cho = self.genessorted[0].split("^")[0]
        pos0 = self.gene2pos[self.genessorted[0]] + self.start_list[self.labels.index(cho)]
        chain = self.gene2chain[self.genessorted[0]]
        # print('pos0,total_size:',pos0, total_size)
        start = self.to_radian(pos0, total_size)
        self.plot_arc_block(start, chain, self.radius_a - shiftlevel * self.shiftratio * self.blockthick)
        gene2location[self.genessorted[0]]= start
        gene2shift[self.genessorted[0]] = shiftlevel
        laststart = start

        for i in range(len(self.genessorted)):
            chi = self.genessorted[i].split("^")[0]
            posi = self.gene2pos[self.genessorted[i]] + self.start_list[self.labels.index(chi)]
            chain = self.gene2chain[self.genessorted[i]]
            # print('pos0,total_size:',pos0, total_size)# pos0,total_size: 2120381 466090045.0
            start = self.to_radian(posi, total_size)
            if(start-laststart < self.blocklength):
                shiftlevel = shiftlevel + 1
            else:
                 shiftlevel=0

            shiftradius = self.radius_a - shiftlevel * self.shiftratio * self.blockthick

        #    print "draw block: ", start, chain, shiftlevel, self.shiftratio, self.blockthick, shiftradius, self.radius_a

        #    self.plot_arc_block(start, chain, self.radius_a - (shiftlevel+1) * self.shiftratio * self.blockthick, "blue")
            gene2location[self.genessorted[i]] = start
            gene2shift[self.genessorted[i]] = shiftlevel
            # print("my start+++++++++++++", self.genessorted[i], chi, gene2pos[self.genessorted[i]], shiftlevel)

            laststart = start
        ### draw Ks lines between genes
        # print("***************************************************")
        for genepair in self.genepair2Ks.keys():
        #    print genepair
            gene1, gene2 = genepair.split(' ')
            Ks = self.genepair2Ks[genepair]
            if (gene1 not in gene2location.keys()):
                # print(gene1)
                continue
            start1 = gene2location[gene1]
            # print(start1)
            shift1 = gene2shift[gene1]
            if (gene2 not in gene2location.keys()):
                # print(gene2)
                continue
            start2 = gene2location[gene2]
            # print(start2)
            shift2 = gene2shift[gene2]
            # print ("start", start1, start2)# start 6.015783710157899 6.133748376342092
            # python abs 绝对值
            # if (Ks < 0.5):
            #     if (abs(start1-start2) < 0.05):
            #         makers = '-'
            #         color = "gray"
            #         self.plot_bez_Ks2(start1, self.radius_a - (shift1+1) * self.shiftratio * self.blockthick, start2, self.radius_a - (shift2+1) * self.shiftratio * self.blockthick, color, makers)
            color, makers, class1 = self.run_plot(Ks, start1, start2)
            if (color != None):
                # print(a0,b0)
                self.plot_bez_Ks2(start1, self.radius_a - (shift1+1) * self.shiftratio * self.blockthick, start2, self.radius_a - (shift2+1) * self.shiftratio * self.blockthick, color, makers)
                lt = [gene1, gene2, start1, shift1, start2, shift2, class1]
                with open(self.savecsv, "a", newline='', encoding='utf-8') as file:
                    writer = csv.writer(file ,delimiter=',')
                    writer.writerow(lt)

            else :
                # print('间距不达标')
                pass
            # # abs(start1-start2) < 0.5 控制连线的长度
            # # Ks控制连线数量
            # if(Ks >=1 and Ks < 1.5 ):
            #      color = "blue"
            #      makers = '-'
            #      if abs(start1-start2) < 0.005:
            #         self.plot_bez_Ks2(start1, self.radius_a - (shift1+1) * self.shiftratio * self.blockthick, start2, self.radius_a - (shift2+1) * self.shiftratio * self.blockthick, color, makers)
            #      #else:
            #     #    plot_bez_Ks(start1, self.radius_a - (shift1+1) * self.shiftratio * self.blockthick, start2, self.radius_a - (shift2+1) * self.shiftratio * self.blockthick, color, makers)
            # if(Ks >= 0.2 and Ks < 1 ):
            #      color = "green"
            #      makers = '-'
            #      if abs(start1-start2) < 0.05:
            #         print(self.radius_a - (shift1+1) * self.shiftratio * self.blockthick,'***', self.radius_a - (shift2+1) * self.shiftratio * self.blockthick)
            #         self.plot_bez_Ks2(start1, self.radius_a - (shift1+1) * self.shiftratio * self.blockthick, start2, self.radius_a - (shift2+1) * self.shiftratio * self.blockthick, color, makers)
            # if(Ks >=0.1 and Ks < 0.11 ):
            #      color = "blue"
            #      makers = '-'
            #      makers1 = '.'
            #      if abs(start1-start2) < 0.05:
            #         self.plot_bez_Ks2(start1, self.radius_a - (shift1+1) * self.shiftratio * self.blockthick, start2, self.radius_a - (shift2+1) * self.shiftratio * self.blockthick, color, makers)
            #      else:
            #         # print('1')
            #         self.plot_bez_Ks(start1, self.radius_a - (shift1+1) * self.shiftratio * self.blockthick, start2, self.radius_a - (shift2+1) * self.shiftratio * self.blockthick, color, makers1)
            # if( Ks < 0.1):
            #      color = "red"
            #      makers = '-'
            #      if abs(start1-start2) < 0.002:
            #         self.plot_bez_Ks3(start1, self.radius_a - (shift1+1) * self.shiftratio * self.blockthick, start2, self.radius_a - (shift2+1) * self.shiftratio * self.blockthick, color, makers)
                 

                 #else:
                    #plot_bez_Ks(start1, self.radius_a - (shift1+1) * self.shiftratio * self.blockthick, start2, self.radius_a - (shift2+1) * self.shiftratio * self.blockthick, color, makers)
                    #    elif(Ks < 0.3):
        #         color = "green"
        #         makers = '-'
        #         if abs(start1-start2) < 0.5:
        #            plot_bez_Ks2(start1, self.radius_a - (shift1+1) * self.shiftratio * self.blockthick, start2, self.radius_a - (shift2+1) * self.shiftratio * self.blockthick, color, makers)
        #         else:
        #            plot_bez_Ks(start1, self.radius_a - (shift1+1) * self.shiftratio * self.blockthick, start2, self.radius_a - (shift2+1) * self.shiftratio * self.blockthick, color, makers)
        #    else:
        #         color = "blue"
        #         makers = '-'
                 #plot_bez_Ks(start1, self.radius_a - (shift1+1) * self.shiftratio * self.blockthick, start2, self.radius_a - (shift2+1) * self.shiftratio * self.blockthick, color, makers)

        ### Again: define shift level and draw gene blocks
        shiftlevel = 0
        gene2shift = {}
        gene2location = {}
        cho = self.genessorted[0].split("^")[0]
        pos0 = self.gene2pos[self.genessorted[0]] + self.start_list[self.labels.index(cho)]
        chain = self.gene2chain[self.genessorted[0]]
        start = self.to_radian(pos0, total_size)
        self.plot_arc_block(start, chain, self.radius_a - shiftlevel * self.shiftratio * self.blockthick)
        #gene2location[self.genessorted[0]]= start
        #gene2shift[self.genessorted[0]] = shiftlevel
        laststart = start

        for i in range(len(self.genessorted)):
            chi = self.genessorted[i].split("^")[0]
            posi = self.gene2pos[self.genessorted[i]] + self.start_list[self.labels.index(chi)]
            chain = self.gene2chain[self.genessorted[i]]
            start = self.to_radian(posi, total_size)
            if(start-laststart < self.blocklength):
                shiftlevel = shiftlevel + 1
            else:
                 shiftlevel=0

            shiftradius = self.radius_a - shiftlevel * self.shiftratio * self.blockthick
        #    print "draw block: ", start, chain, shiftlevel, self.shiftratio, self.blockthick, shiftradius, self.radius_a
            self.plot_arc_block(start, chain, self.radius_a - (shiftlevel+1) * self.shiftratio * self.blockthick)
        #    gene2location[self.genessorted[i]] = start
        #    gene2shift[self.genessorted[i]] = shiftlevel
        #    print self.genessorted[i], start

            laststart = start


        #### draw tick showing scale of chromosomes
        for i in range(chr_number):
           pos = 0
           while pos < self.chro2len[self.labels[i]]:
              xx1, yy1 = self.transform_pt(i, int(pos), self.radius_b, total_size)
              xx2, yy2 = self.transform_pt(i, int(pos), self.radius_b-0.003, total_size)
              plot([xx1, xx2], [yy1, yy2], "k-", lw = .2)
              pos = pos + 10000000

        #
        root.set_xlim(-.5, .5)
        root.set_ylim(-.5, .5)
        root.set_axis_off()
        savefig(self.savefile, dpi=500)
        sys.exit(0)
