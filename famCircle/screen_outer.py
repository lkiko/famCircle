###
# -*- coding: UTF-8 -*-

import re
import os
from math import *
import csv
import pandas as pd
from matplotlib.patches import *
from pylab import *
from famCircle.bez import *


class screen_outer():
    def __init__(self, options):
        for k, v in base_conf:
            setattr(self, str(k), v)
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)

