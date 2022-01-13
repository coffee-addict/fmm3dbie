import os
import sys
import warnings
import re
import time
from datetime import datetime
from glob import glob
import numpy as np
import matplotlib.pyplot as plt

f_input = './geometry.txt'

def main():
    f = open(f_input, mode='r', encoding='utf-8')
    n = 0
    for i, line in enumerate(f):
        n += 1
    f.close()
    x = [0]*n
    y = [0]*n
    z = [0]*n
    f = open(f_input, mode='r', encoding='utf-8')
    for i, line in enumerate(f):
        columns = line.strip().strip('\r').strip('\n').split('  ')
        columns = list(map(lambda s: float(s.strip()), columns))
        x[i],y[i],z[i] = columns
    f.close()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z)
#    ax.plot_trisurf(x, y, z, linewidth=0.2,antialiased=True)
#    ax.plot_surface(x, y, z)

    lim_axes_min = min(min(x),min(y),min(z))
    lim_axes_max = max(max(x),max(y),max(z))

    ax.set_xticks(np.linspace(lim_axes_min,lim_axes_max,5))
    ax.set_yticks(np.linspace(lim_axes_min,lim_axes_max,5))
    ax.set_zticks(np.linspace(lim_axes_min,lim_axes_max,5))
    ax.axes.set_xlim3d(left=lim_axes_min, right=lim_axes_max)
    ax.axes.set_ylim3d(bottom=lim_axes_min, top=lim_axes_max)
    ax.axes.set_zlim3d(bottom=lim_axes_min, top=lim_axes_max)


    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()
    
if __name__ == '__main__': main()
