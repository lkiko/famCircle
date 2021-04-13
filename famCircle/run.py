# -*- coding: UTF-8 -*-
import argparse
import os
import sys
import configparser
import pandas as pd
import famCircle
import famCircle.bez as bez
from famCircle.hmmer import hmmer
from famCircle.lookKs import lookKs
from famCircle.circle import circle
from famCircle.outer import outer
from famCircle.inner import inner

parser = argparse.ArgumentParser(
    prog = 'famCircle', usage = '%(prog)s [options]', epilog = "", formatter_class = argparse.RawDescriptionHelpFormatter,)
parser.description = '''\
圈图构建
    -------------------------------------- '''
parser.add_argument("-v", "--version", action = 'version', version='0.1.1')
parser.add_argument("-hmmer", dest = "hmmer",
                    help = "基因家族鉴定")
parser.add_argument("-lookKs", dest = "lookKs",
                    help = "Ks可视化")
parser.add_argument("-circle", dest = "circle",
                    help = "共线性可视化")
parser.add_argument("-outer", dest = "outer",
                    help = "当基因组复杂且重复基因多时在圈图外围显示，构建发射状基因圈图")
parser.add_argument("-inner", dest = "inner",
                    help = "当基因组简单且重复基因比较少时在圈图内显示，构建内卷状基因圈图")

args = parser.parse_args()

def run_hmmer():
    options = bez.load_conf(args.hmmer, 'hmmer')
    hmmer1 = hmmer(options)
    hmmer1.run()

def run_lookKs():
    options = bez.load_conf(args.lookKs, 'lookKs')
    lookKs1 = lookKs(options)
    lookKs1.run()

def run_circle():
    options = bez.load_conf(args.circle, 'circle')
    circle1 = circle(options)
    circle1.run()

def run_outer():
    options = bez.load_conf(args.outer, 'outer')
    outer1 = outer(options)
    outer1.run()

def run_inner():
    options = bez.load_conf(args.inner, 'inner')
    inner1 = inner(options)
    inner1.run()

def module_to_run(argument):
    switcher = {
        'hmmer': run_hmmer,
        'lookKs': run_lookKs,
        'circle': run_circle,
        'outer': run_outer,
        'inner': run_inner,
    }
    return switcher.get(argument)()


def main():
    path = famCircle.__path__[0]
    options = {
               'hmmer': 'hmmer.conf',
               'lookKs': 'lookKs.conf',
               'circle': 'circle.conf',
               'outer': 'outer.conf',
               'inner': 'inner.conf',
               }
    for arg in vars(args):
        value = getattr(args, arg)
        # print(value)
        if value is not None:
            if value in ['?', 'help', 'example']:
                f = open(os.path.join(path, 'example', options[arg]))
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

