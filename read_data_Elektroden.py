__author__ = 'plappert'
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import glob

# Code, um in csv umgewandeltes ods sheet zu lesen. Das waren die Daten, wo ich ganz am Anfang per Hand das Feld der beiden Elektroden (E1, E2) ausgemessen hab

def read_data(filename): #Funktion liest die einzelnen Zeilen des csv sheets und macht aus strings Zahlen und loescht die Kommas zwischen den Zahlen
    f = file(filename)
    lines = f.readlines()
    f.close()
    XPos = []
    YPos = []
    amps = None
    line_counter = 0
    for l in lines:
        if 'Elektrode' in l:
            continue
        if 'Position' in l:
            XPos = map(float, l.strip().split(',')[1:])
        else:
            temp = map(float, l.strip().split(','))
            YPos.append(temp[0])
            if amps is None:
                amps = np.zeros((len(lines)-2, len(temp)-1))
            amps[line_counter, :] = temp[1:]
            line_counter += 1
    return XPos, YPos, amps


def heat_map(x_pos, y_pos, amplitude): # Fuktion erstellt fuer jede Elektrode eine Heatmap über den Befehl imshow, über extent koennen die x und die y achse bestimmt werden
   plt.imshow(amplitude, interpolation='bicubic', cmap=cm.jet, vmin=.0, vmax=np.max(amplitude), extent=[0,105,0,55])
   plt.xlabel('X Position [cm]')
   plt.ylabel('Y Position [cm]')
   cb = plt.colorbar() # erstellt eine Legende zu den verschiedenen Farbintensitaeten
   cb.set_label('mV') # beschriftung der legende
   if 'elektrode1' in fl:
       plt.title('Feldintensitaet Elektrode 1')
       plt.savefig('2015-06-23_Heatmap_Elektrode1.pdf')
   if 'elektrode2' in fl:
       plt.title('Feldintensitaet Elektrode 2')
       plt.savefig('2015-06-23_Heatmap_Elektrode2.pdf')
   plt.show()


if __name__ == '__main__':
    files = glob.glob('/home/plappert/Documents/2015-06-23elektrode*') #die beider csv sheets (pro elektrode eins) werden geoeffnet
    print files
    for fl in files:
        x_pos, y_pos, amplitude = read_data(fl)
        heat_map(x_pos, y_pos, amplitude)






