# -*- encoding: utf-8 -*-
'''
@File        :run.py
@Time        :2021/09/28 09:04:51
@Author        :charles kiko
@Version        :1.0
@Contact        :charles_kiko@163.com
@Desc        :主程序
'''


import argparse
import os
import sys
import configparser
import pandas as pd
import famCircle
import famCircle.bez as bez
from famCircle.filterWGD import filterWGD
from famCircle.hmmer import hmmer
from famCircle.screen import screen
from famCircle.Ks import Ks
from famCircle.Ks_block import Ks_block
from famCircle.circle import circle
from famCircle.line import line
from famCircle.sline import sline
from famCircle.lineblock import lineblock
from famCircle.circle_all import circle_all
from famCircle.loops import loops
from famCircle.re_loops import re_loops
from famCircle.loops_ks import loops_ks
from famCircle.circle_family import circle_family
from famCircle.famloc import famloc
from famCircle.filter_unit import filter_unit
from famCircle.abs_tandem import abs_tandem

from famCircle.outer import outer
from famCircle.inner import inner
from famCircle.part_out import part_out
from famCircle.cdotplot import cdotplot
from famCircle.typing import typing
from famCircle.part import part
from famCircle.family_pair import family_pair

from famCircle.classification import classification

from famCircle.trd import trd

from famCircle.msa import msa
from famCircle.tree import tree

parser = argparse.ArgumentParser(
    prog = 'famCircle', usage = '%(prog)s [options]', epilog = "", formatter_class = argparse.RawDescriptionHelpFormatter,)
parser.description = '''\
runing famCircle
    -------------------------------------- '''
parser.add_argument("-v", "--version", action = 'version', version='0.2.5')

parser.add_argument("-fu", dest = "filter_unit",
                    help = "Filter blast files to keep the best match that is not self.")
parser.add_argument("-ca", dest = "circle_all",
                    help = "Draws a collinearity circle plot, with support for WGDIDraws a collinearity circle plot, with support for MCScanX/famCircle/CollinearScan/WGDI file formats.")
parser.add_argument("-lp", dest = "loops",
                    help = "The genomic collinearity circle map was drawn and the GC content of repetitive sequences was analyzed.")
parser.add_argument("-lpks", dest = "loops_ks",
                    help = "The genomic collinearity circle map with KS was drawn and the content of GC in repeat sequence was analyzed.")
parser.add_argument("-ilp", dest = "re_loops",
                    help = "The genomic collinearity circle map was drawn and the GC content of repetitive sequences was analyzed.")
parser.add_argument("-l", dest = "line",
                    help = "Draw parallel collinear diagrams.")
parser.add_argument("-sl", dest = "sline",
                    help = "Draw parallel collinear diagrams.")
parser.add_argument("-lb", dest = "lineblock",
                    help = "To plot the linear stack of genes.")
parser.add_argument("-cd", dest = "cdotplot",
                    help = "Plots the genome collinearity matrix, which only supports collinear files.")
parser.add_argument("-c", dest = "circle",
                    help = "Draw a collinearity circle plot with ks.")
parser.add_argument("-ks", dest = "Ks",
                    help = "Calculate ka/ks, multi-process acceleration version.")
parser.add_argument("-kb", dest = "Ks_block",
                    help = "Ks kernel density plot of genomic collinear blocks.")
parser.add_argument("-hmm", dest = "hmmer",
                    help = "Scanning gene family members with HMM model.")
parser.add_argument("-fp", dest = "family_pair",
                    help = "Extract alignment records about family members from BLAST/MCScanX files.")
parser.add_argument("-cf", dest = "circle_family",
                    help = "Chromosomal localization map of gene family members.")
parser.add_argument("-fl", dest = "famloc",
                    help = "Chromosomal localization map of gene family members.")

parser.add_argument("-s", dest = "screen",
                    help = "To map the chromosome locations of gene families.")
parser.add_argument("-at", dest = "abs_tandem",
                    help = "Absolute tandem repeat identification of genome.")
parser.add_argument("-o", dest = "outer",
                    help = "Gene family circle diagram, external accumulation diagram.")
parser.add_argument("-po", dest = "part_out",
                    help = "Chromosomal Localization of Ring Gene Family Stack.")
parser.add_argument("-i", dest = "inner",
                    help = "Gene family circle graph, internal accumulation graph.")
parser.add_argument("-p", dest = "part",
                    help = "The linear gene family accumulation map on the chromosome.")
parser.add_argument("-msa", dest = "msa",
                    help = "Alignment of multiple sequences within a gene cluster.")
parser.add_argument("-tree", dest = "tree",
                    help = "Constructing a phylogenetic tree for the gene cluster using software PhyML/IQtree/fasttree.")


parser.add_argument("-t", dest = "typing",
                    help = "Temporarily unavailable.")
parser.add_argument("-f", dest = "filterWGD",
                    help = "Temporarily unavailable.")
parser.add_argument("-class", dest = "classification",
                    help = "Temporarily unavailable.")
parser.add_argument("-trd", dest = "trd",
                    help = "Temporarily unavailable.")



args = parser.parse_args()

def run_filter_unit():
    options = bez.load_conf(args.filter_unit, 'filter_unit')
    filter_unit1 = filter_unit(options)
    filter_unit1.run()

def run_filterWGD():
    options = bez.load_conf(args.filterWGD, 'filterWGD')
    filterWGD1 = filterWGD(options)
    filterWGD1.run()

def run_classification():
    options = bez.load_conf(args.classification, 'classification')
    classification1 = classification(options)
    classification1.run()

def run_trd():
    options = bez.load_conf(args.trd, 'trd')
    trd1 = trd(options)
    trd1.run()

def run_cdotplot():
    options = bez.load_conf(args.cdotplot, 'cdotplot')
    cdotplot1 = cdotplot(options)
    cdotplot1.run()

def run_hmmer():
    options = bez.load_conf(args.hmmer, 'hmmer')
    hmmer1 = hmmer(options)
    hmmer1.run()

def run_screen():
    options = bez.load_conf(args.screen, 'screen')
    screen1 = screen(options)
    screen1.run()

def run_msa():
    options = bez.load_conf(args.msa, 'msa')
    msa1 = msa(options)
    msa1.run()

def run_tree():
    options = bez.load_conf(args.tree, 'tree')
    tree1 = tree(options)
    tree1.run()

def run_Ks():
    options = bez.load_conf(args.Ks, 'Ks')
    lookKs0 = Ks(options)
    lookKs0.run()

def run_Ks_block():
    options = bez.load_conf(args.Ks_block, 'Ks_block')
    lookKs0 = Ks_block(options)
    lookKs0.run()

def run_abs_tandem():
    options = bez.load_conf(args.abs_tandem, 'abs_tandem')
    abs_tandem0 = abs_tandem(options)
    abs_tandem0.run()

def run_typing():
    options = bez.load_conf(args.typing, 'typing')
    typing1 = typing(options)
    typing1.run()

def run_circle():
    options = bez.load_conf(args.circle, 'circle')
    circle1 = circle(options)
    circle1.run()

def run_loops():
    options = bez.load_conf(args.loops, 'loops')
    loops1 = loops(options)
    loops1.run()

def run_loops_ks():
    options = bez.load_conf(args.loops_ks, 'loops_ks')
    loops_ks1 = loops_ks(options)
    loops_ks1.run()

def run_re_loops():
    options = bez.load_conf(args.re_loops, 're_loops')
    re_loops1 = re_loops(options)
    re_loops1.run()

def run_line():
    options = bez.load_conf(args.line, 'line')
    circle1 = line(options)
    circle1.run()

def run_sline():
    options = bez.load_conf(args.sline, 'sline')
    circle1 = sline(options)
    circle1.run()

def run_lineblock():
    options = bez.load_conf(args.lineblock, 'lineblock')
    lineblock1 = lineblock(options)
    lineblock1.run()

def run_circle_all():
    options = bez.load_conf(args.circle_all, 'circle_all')
    circle0 = circle_all(options)
    circle0.run()

def run_circle_family():
    options = bez.load_conf(args.circle_family, 'circle_family')
    circle0 = circle_family(options)
    circle0.run()

def run_famloc():
    options = bez.load_conf(args.famloc, 'famloc')
    circle0 = famloc(options)
    circle0.run()

def run_outer():
    options = bez.load_conf(args.outer, 'outer')
    outer1 = outer(options)
    outer1.run()

def run_part_out():
    options = bez.load_conf(args.part_out, 'part_out')
    part_out1 = part_out(options)
    part_out1.run()

def run_inner():
    options = bez.load_conf(args.inner, 'inner')
    inner1 = inner(options)
    inner1.run()

def run_part():
    options = bez.load_conf(args.part, 'part')
    inner1 = part(options)
    inner1.run()

def run_family_pair():
    options = bez.load_conf(args.family_pair, 'family_pair')
    family_pair1 = family_pair(options)
    family_pair1.run()

def module_to_run(argument):
    switcher = {
        'filter_unit': run_filter_unit,
        'filterWGD': run_filterWGD,
        'classification': run_classification,
        'cdotplot': run_cdotplot,

        'trd': run_trd,

        'hmmer': run_hmmer,
        'abs_tandem': run_abs_tandem,
        'screen': run_screen,
        'Ks': run_Ks,
        'Ks_block': run_Ks_block,
        'typing': run_typing,
        'circle': run_circle,
        'line': run_line,
        'sline': run_sline,
        'lineblock': run_lineblock,
        'circle_all': run_circle_all,
        'loops': run_loops,
        're_loops': run_re_loops,
        'loops_ks': run_loops_ks,
        'circle_family': run_circle_family,
        'famloc': run_famloc,

        'outer': run_outer,
        'part_out': run_part_out,
        'inner': run_inner,
        'part': run_part,
        'msa': run_msa,
        'tree': run_tree,
        'family_pair': run_family_pair,
    }
    return switcher.get(argument)()

def main():
    path = famCircle.__path__[0]
    options = {
               'filter_unit': 'filter_unit.conf',
               'filterWGD': 'filterWGD.conf',
               'classification': 'classification.conf',
               'trd': 'trd.conf',
               'cdotplot': 'cdotplot.conf',
               'hmmer': 'hmmer.conf',
               'screen': 'screen.conf',
               'Ks': 'Ks.conf',
               'Ks_block': 'Ks_block.conf',
               'typing': 'typing.conf',
               'circle': 'circle.conf',
               'line': 'line.conf',
               'sline': 'sline.conf',
               'lineblock': 'lineblock.conf',
               'circle_all': 'circle_all.conf',
               'loops': 'loops.conf',
               're_loops': 're_loops.conf',
               'loops_ks': 'loops_ks.conf',
               'circle_family': 'circle_family.conf',
               'famloc': 'famloc.conf',

               'outer': 'outer.conf',
               'part_out': 'part_out.conf',
               'abs_tandem': 'abs_tandem.conf',
               'inner': 'inner.conf',
               'part': 'part.conf',
               'msa': 'msa.conf',
               'tree': 'tree.conf',
               'family_pair': 'family_pair.conf',
               }
    for arg in vars(args):
        value = getattr(args, arg)
        # print(value)
        if value is not None:
            if value in ['?', 'help', 'example']:
                f = open(os.path.join(path, 'example', options[arg]), encoding='utf-8')
                print(f.read())
            elif value == 'e':
                out = '''\
        File example
        [fpchrolen]
        chromosomes number_of_bases
        *   *
        *   *
        *   *
        [fpgff]
        chromosomes gene    start   end
        *   *   *   *
        *   *   *   *
        *   *   *   *
        [fpgenefamilyinf]
        gene1   gene2   Ka  Ks
        *   *   *   *
        *   *   *   *
        *   *   *   *
        [alphagenepairs]
        gene1   gene2
        *   *   *
        *   *   *
        *   *   *

        The file columns are separated by Tab
        -----------------------------------------------------------    '''
                print(out)
            elif not os.path.exists(value):
                print(value+' not exits')
                sys.exit(0)
            else:
                module_to_run(arg)

