__author__ = 'plappert'
__author__ = 'plappert'
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from IPython import embed
import glob
import nix
from scipy.stats import linregress
import itertools
from pylab import *
import math


fig, ax1 = plt.subplots()

fs = 44100
t = np.arange(-0.002, .002, 1.0/fs)
f0 = 900
phi = np.pi/2
A = 1
x = A * np.sin(2 * np.pi * f0 * t + phi)

ax1.plot(t, x, linewidth=3)
plt.axis([-0.002, .002, -1, 1.2])
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax1.get_yticklabels(), visible=False)
plt.savefig('Frequenz1.pdf')
plt.show()

fig, ax2 = plt.subplots()
fs = 44100
t = np.arange(-0.002, .002, 1.0/fs)
f1 = 450
phi = np.pi/2
A = 1
x = A * np.sin(2 * np.pi * f1 * t + phi)

ax2.plot(t, x, linewidth=3)
plt.axis([-0.002, .002, -1, 1.2])
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
plt.savefig('Frequenz2.pdf')
plt.show()