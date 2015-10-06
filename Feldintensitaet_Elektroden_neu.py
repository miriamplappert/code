__author__ = 'plappert'
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import glob
import pyrelacs
import gc
import scipy.io as scio

# die via relacs gesammelten Daten, um die felder der beiden elektroden (E1 und E2) auszumessen werden hier verarbeitet und es wird fuer jede Elektrode eine heatmap generiert

def merge_amplitudes(matrix, pos_amp, dx, dy):
    for p in pos_amp:
        x_index = p[0][0] / dx - 1
        y_index = p[0][1] / dy - 1
        if matrix[y_index, x_index] == 0.0:
            matrix[y_index, x_index] = p[1]
    return matrix


def convert_data(data):
    time = np.zeros(len(data))
    voltage = np.zeros(len(data))
    for index,d in enumerate(data):
        time[index], voltage[index] = map(float, d.strip().split())
    return time, voltage


def load_data(datei):
    print datei
    data = pyrelacs.DataClasses.load(datei)
    m,k,d = data.select({('output channel',): '0LocalEField-1'})
    electrode_1 = analyse_data(m,d)
    m,k,d = data.select({('output channel',): '0LocalEField-2'})
    electrode_2 = analyse_data(m,d)
    return electrode_1, electrode_2


def analyse_data(meta, data):
    pos_amp = []
    for index in range(len(meta)):
        m = meta[index]
        d = data[index]
        time, voltage = convert_data(d)
        mittlere_amplitude = np.percentile(voltage, 90) - np.percentile(voltage, 10)
        position = (float(m['x-pos']), float (m['y-pos']))
        pos_amp.append((position,mittlere_amplitude))
    return pos_amp

def heat_map(amplitude, title):

    plt.imshow((amplitude), interpolation='bicubic', cmap=cm.jet, vmin=.01, vmax=np.max(amplitude), extent=[0,105,0,55])
    plt.xlabel('X Position [cm]')
    plt.ylabel('Y Position [cm]')
    cb = plt.colorbar()
    cb.set_label('Amplitude [mV]')
    plt.title(title)
    plt.savefig('Heatmap_' + title + '.pdf')
    plt.close()

def contourmap(amplitude, title):
    cs = plt.contour(amplitude)
    plt.xlabel('X Position [cm]')
    plt.ylabel('Y Position [cm]')
    cb = plt.colorbar()
    cb.set_label('Amplitude [mV]')
    plt.clabel(cs,inline=1,fontsize=10)
    plt.title(title)
    plt.savefig('Contourmap_' + title + '.pdf')
    plt.show()


if __name__ == '__main__':
    amplituden_elektrode1 = np.zeros((9,19))
    amplituden_elektrode2 = np.zeros((9,19))
    saved_data = glob.glob('heatmap_data.mat')
    if len(saved_data) > 0:
        heatmap_data = scio.loadmat('heatmap_data.mat')
        amplituden_elektrode1 = heatmap_data['matrix_electrode_1']
        amplituden_elektrode2 = heatmap_data['matrix_electrode_2']
    else:
        files = glob.glob('/home/plappert/Data/2015-06-30*')
        dx = 5
        dy = 5
        for f in files:
            dateien = glob.glob(f+'/localField*.dat')
            for d in dateien:
                electrode_1, electrode_2 = load_data(d)
                amplituden_elektrode1 = merge_amplitudes(amplituden_elektrode1,electrode_1, dx, dy)
                amplituden_elektrode2 = merge_amplitudes(amplituden_elektrode2,electrode_2, dx, dy)
                gc.collect()
        scio.savemat('heatmap_data.mat', {'electrode_1': electrode_1, 'electrode_2':electrode_2,
                                          'matrix_electrode_1': amplituden_elektrode1,
                                          'matrix_electrode_2': amplituden_elektrode2})
    heat_map(amplituden_elektrode1, 'electrode 1')
    heat_map(amplituden_elektrode2, 'electrode 2')

    contourmap(amplituden_elektrode1, 'electrode 1')
    contourmap(amplituden_elektrode2, 'electrode 2')
