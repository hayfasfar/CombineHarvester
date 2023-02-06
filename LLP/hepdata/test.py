import matplotlib
matplotlib.use('Agg')
import style
import ROOT
ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetPalette(ROOT.kDarkRainBow)
import json
from array import array
from scipy import interpolate
import numpy as np
import pandas as pd
import os


hist_displaced_dirac = get_graph("hepdata/HEPData-ins1736526-v1.root", "Table 3/Graph1D_y1")

