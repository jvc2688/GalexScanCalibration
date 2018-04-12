import c3
from astropy.io import fits as pyfits
import numpy as np
from astropy import wcs as pywcs
import matplotlib.pyplot as plt
import sys
import aplpy
import os
from sklearn.neighbors import KernelDensity
import csv
import math

fig = aplpy.FITSFigure('../fits/count_map_05_gPr_res_r_flip.fits')
fig.show_grayscale(invert=True)
fig.save('../plots/corr_10/test.png')