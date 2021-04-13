# -*- coding: UTF-8 -*-
import configparser
import os
import re
import famCircle
import numpy as np
import pandas as pd
from Bio import Seq, SeqIO, SeqRecord
import codecs
# Bezier functions

def config():
    conf = configparser.ConfigParser()
    conf.read(os.path.join(famCircle.__path__[0], 'conf.ini'))
    return conf.items('ini')

def load_conf(file, section):
    conf = configparser.ConfigParser()
    conf.read(file)
    return conf.items(section)

def calculate_coef(p0, p1, p2, p3):
    c = 3*(p1 - p0)
    b = 3*(p2 - p1) -c
    a = p3 - p0 - c - b
    return c, b, a

def Bezier(plist, t):
    # p0 : origin, p1, p2 :control, p3: destination
    p0, p1, p2, p3 = plist
    # calculates the coefficient values
    c, b, a = calculate_coef(p0, p1, p2, p3)
    tsquared = t**2
    tcubic = tsquared*t
    return a*tcubic + b*tsquared + c*t + p0