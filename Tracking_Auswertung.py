__author__ = 'plappert'
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from IPython import embed
import glob
import nix
from read_data_versuch4 import videofiles1, videofiles2, videofiles3, videofiles4, videofiles5, videofiles6
from scipy.stats import linregress
import itertools
from pylab import *

def orientation_to_electrode(orientations, x_positions, y_positions, distance_to_E1, distance_to_E2, k):

    print x_positions

    index_small_distance_E1 = []
    for i in np.arange(len(distance_to_E1)): # for schleife erstellt liste mit den indices, an denen im eod array eine 0 steht
        if distance_to_E1[i] < (100*0.12):
            index_small_distance_E1.append(i)

    index_small_distance_E2 = []
    for i2 in np.arange(len(distance_to_E2)): # for schleife erstellt liste mit den indices, an denen im eod array eine 0 steht
        if distance_to_E2[i2] < (100*0.12):
            index_small_distance_E2.append(i2)


    x_position_near_e1 = [i for j, i in enumerate(x_positions) if j not in index_small_distance_E1]
    y_position_near_e1 = [i for j, i in enumerate(y_positions) if j not in index_small_distance_E1]
    orientation_near_e1 = [i for j, i in enumerate(orientations) if j not in index_small_distance_E1]

    x_position_near_e2 = [i2 for j2, i2 in enumerate(x_positions) if j2 not in index_small_distance_E2]
    y_position_near_e2 = [i2 for j2, i2 in enumerate(y_positions) if j2 not in index_small_distance_E2]
    orientation_near_e2 = [i2 for j2, i2 in enumerate(orientations) if j2 not in index_small_distance_E2]



    return

def distance_velocity_plot(E1_distance, E2_distance, velocity, filename):

    m,y_achsenabschnitt, r_value, p_value, std_err = linregress(E1_distance[:-1], velocity)
    plt.plot(E1_distance[:-1], velocity, 'yo', label='Pearson\'r = %.3f, p = %.2e'%(r_value, p_value))
    y =[] #the following lines create the regression line
    for t in E1_distance[:-1]:
        f = m * t + y_achsenabschnitt
        y.append(f)
    plt.plot(E1_distance[:-1], y) #plots regression
    plt.title('Geschwindigkeit in Abhaengigkeit zum Elektrodenabstand ' + filename)
    plt.xlabel('Distanz [Zentimeter]')
    plt.ylabel('Geschwindigkeit [Zentimeter pro Sekunde]')
    plt.xlim(0,80)
    plt.ylim(0,70)
    plt.legend(loc=2, numpoints=1, markerscale=0., frameon=False)
    plt.savefig('distance_velocity_plot' + filename + 'E1.pdf')
    plt.close()

    m2,y_achsenabschnitt2, r_value2, p_value2, std_err2 = linregress(E2_distance[:-1], velocity)
    plt.plot(E2_distance[:-1], velocity, 'yo', label='Pearson\'r = %.3f, p = %.2e'%(r_value2, p_value2))
    y2 =[] #the following lines create the regression line
    for t2 in E2_distance[:-1]:
        f2 = m2 * t2 + y_achsenabschnitt2
        y2.append(f2)
    plt.plot(E2_distance[:-1], y2) #plots regression
    plt.title('Geschwindigkeit in Abhaengigkeit zum Elektrodenabstand ' + filename)
    plt.xlabel('Distanz [Pixel]')
    plt.ylabel('Geschwindigkeitx[Pixel pro Sekunde]')
    plt.xlim(0,80)
    plt.ylim(0,70)
    plt.legend(loc=2, numpoints=1, markerscale=0., frameon=False)
    plt.savefig('distance_velocity_plot' + filename + 'E2.pdf')
    plt.close()


def small_distance_amount_bar_plot(small_E1_distance_amount, small_E2_distance_amount, filename):
    '''
    function creates a bar plot that shows how often the fish was nearer than 100 Pixel to electrode1 and 2
    :param small_E1_distance_amount: amount of measurings where distance between fish and electrode1 was smaller than 100 pixel
    :param small_E2_distance_amount: amount of measurings where distance between fish and electrode2 was smaller than 100 pixel
    :param filename:
    :return:bar plot
    '''
    N = 1
    ind = np.arange(N)  # the x locations for the groups
    width = 0.1  # the width of the bars
    fig, ax = plt.subplots()

    ## the bars
    first_bar = ax.bar(ind, small_E1_distance_amount, width, color='green')
    second_bar = ax.bar(ind + width, small_E2_distance_amount, width, color='red')


    ax.set_ylabel('Anzahl Naeherungen an Elektrode')
    ax.set_title('Elektrodenannaeherung ' + filename)
    ax.set_xticks(ind+width)
    ax.legend((first_bar[0], second_bar[0]), ('Elektrode 1', 'Elektrode 2'))
    plt.savefig('small_distance_to_electrode_abundance_' + filename + '.pdf')
    plt.close()


def small_distance_to_electrode_abundance(E1_distance, E2_distance):
    '''
    function calculates how often the fish gets nearer than 100 Pixels to electrode1 bzw. 2
    :param E1_distance: list of distances to electrode1
    :param E2_distance: list of distances to electrode2
    :return: two lists that contain how many times the fish was nearer than 100 Pixel to electrode1 bzw. 2
    '''

    small_E1_distance = []
    small_E2_distance = []

    for i in np.arange(len(E1_distance)):

        if E1_distance[i] < (100*0.12):
            small_E1_distance.append(E1_distance[i])

        if E2_distance[i] < (100*0.12):
            small_E2_distance.append(E2_distance[i])

        else:
            continue

    small_E1_distance_amount = len(small_E1_distance)
    small_E2_distance_amount = len(small_E2_distance)

    return small_E1_distance_amount, small_E2_distance_amount


def velocity_time_plot(velocity, time, filename):
    '''
    function creates a plot that shows the velocity for every time of one videofile
    :param velocity: list of velocities
    :param time: list of time ticks
    :param filename: name of file
    :return: velocity plot
    '''


    plt.plot(time[:-1], velocity) #time has to get minus one because the lists are not equal
    plt.title('Geschwindigkeit in Abhaengigkeit zur Zeit ' + filename)
    plt.xlabel('Zeit[s]')
    plt.ylabel('Geschwindigkeitx[Zentimeter pro Sekunde]')
    plt.savefig('velocity_time_plot' + filename + '.pdf')
    plt.close()


def get_velocity(x_pos, y_pos, pos_times):
    '''
    function calculates the velocity (speed) of the fish at every time
    :param x_pos: list of x_positions
    :param y_pos: list of y_positions
    :param pos_times: list of time ticks that fit to the positions
    :return: list of velocities for one videofile
    '''

    velocities = []

    for i in np.arange(len(x_pos)):

        if i == len(x_pos) - 1:
            continue
        else:
            distance = (np.sqrt((x_pos[i] - x_pos[i + 1])**2 + (y_pos[i] - y_pos[i + 1])**2)*0.12) #calculates the distance between two positions using the pythagoras hypothesis, *0.12 converts pixel to cm
            time = pos_times[i + 1] - pos_times[i] #calculates time between the two positions
            v = distance / time #calculates velocity after the general formula v=s/t
            velocities.append(v)
    return velocities


def distance_time_plot(distance_to_E1, distance_to_E2, pos_times, filename):
    '''
    function creates two plots, that means for every electrode one. therefore the distance to electrode 1 bzw. electrode2 is plotted against the time
    :param distance_to_E1: list of distances to electrode1
    :param distance_to_E2:list of distances to electrode2
    :param pos_times: list of times that fits to the distances
    :param filename: name of file
    :return: two plots
    '''
    plt.plot(pos_times, distance_to_E1)
    plt.title("Distanz zu Elektrode1 " + filename)
    plt.xlabel('Zeit [s]')
    plt.ylabel('Distanz [Zentimeter]')
    plt.ylim(0, 80)
    plt.savefig('distance_time_plot_' + filename + '_electrode1.pdf')
    plt.close()

    plt.plot(pos_times, distance_to_E2)
    plt.title("Distanz zu Elektrode2 " + filename)
    plt.xlabel('Zeit [s]')
    plt.ylabel('Distanz [Zentimeter]')
    plt.ylim(0, 80)
    plt.savefig('distance_time_plot_'+ filename + '_electrode2.pdf')
    plt.close()


def get_distance_to_electrode(x_positions, y_positions, E1_coordinates, E2_coordinates):
    '''
    function gets distance to electrode 1 and electrode 2 for every single position of a trackingfile
    :param x_positions: list of x_positions of one video
    :param y_positions: list of y_positons of one video
    :param E1_coordinates: tuple that contains the x and the y position of electrode1
    :param E2_coordinates: tuple that contains the x and the y position of electrode2
    :return: returns two lists: distances to electrode1 and distances to electrode2
    '''
    distance_to_E1 = []
    distance_to_E2 = []

    for i in np.arange(len(x_positions)): #for loop goes throug the indices of x_positions that are the same for y_positions
        distance_to_E1.append(np.sqrt((x_positions[i] - E1_coordinates[0])**2 + (y_positions[i] - E1_coordinates[1])**2)*0.12) #calculates distance between fish position and electrode1 by using the pythagoras hypothesis: a**2 + b**2 = c**2
        distance_to_E2.append(np.sqrt((x_positions[i] - E2_coordinates[0])**2 + (y_positions[i] - E2_coordinates[1])**2)*0.12) # *0.12 converts pixel to cm

    return distance_to_E1, distance_to_E2


def analyse_tracking_data(xpos, ypos, keys, pos_time, orientation, E1_coordinates, E2_coordinates):
    '''
    function gets dictionaries xpos, ypos and pos_time from function read data. function goes into the dictionaries and allows a further analysis of every single videofile.
    This function contains many other functions that analyse the distance to the electrodes, the velocity and do some plots for every single videofile
    :param xpos: dictionary that contains the x_positions under the key of the filename of a videofile
    :param ypos: dicitonary that conatains the y_positions
    :param keys: contains the filenames. which are the key to dictionary xpos, ypos and pos_time
    :param pos_time: dictionary that contains the fitting times to the positions, key are also the filenames
    :return: calls many other functions for further analysis
    '''

    for k in keys: #for loop goes through the filenames
        x_positions = xpos[k] # calls list under the respective filename k
        y_positions = ypos[k]
        pos_times = pos_time[k]
        orientations = orientation[k]

        #further analysis functions:

        distance_to_E1, distance_to_E2 = get_distance_to_electrode(x_positions, y_positions, E1_coordinates, E2_coordinates)
        distance_time_plot(distance_to_E1, distance_to_E2, pos_times, k)
        velocities = get_velocity(x_positions, y_positions, pos_times)
        velocity_time_plot(velocities, pos_times, k)
        small_E1_distance_amount, small_E2_distance_amount = small_distance_to_electrode_abundance(distance_to_E1, distance_to_E2)
        small_distance_amount_bar_plot(small_E1_distance_amount, small_E2_distance_amount, k)
        distance_velocity_plot(distance_to_E1, distance_to_E2, velocities, k)
        orientation_to_electrode(orientations, x_positions, y_positions, distance_to_E1, distance_to_E2, k)

def get_h5_filenames(videofiles):
    '''

    :param videofiles: wird aus dem pythonfile 'read_data_versuch 4 geladen'
    :return: funktion gibt den pfad der h5 files des trackings fuer jeden fisch zurueck
    '''
    h5_file_names = []
    for v in videofiles:
        h5_file = '/home/plappert/Data/video_files/' + v[0:(len(v) - 4)] + '/' + v[0:(len(v) - 3)] + 'h5'
        h5_file_names.append(h5_file)
    return h5_file_names


def read_data(filenames):
    '''
    funktion liest die h5 tracking datei aus und gibt die wichtigsten trackingdaten als listen pro fisch zurueck
    :param filenames:
    :return: funktion gibt position und orientierung des fisches mit den jeweiligen zeitpunkten zurueck
    '''


    estimated_xpos = dict()
    estimated_ypos = dict()
    estimated_pos_times = dict()
    estimated_orientations = dict()
    xpos = dict()
    ypos = dict()
    pos_times = dict()
    orientations = dict()
    keys = []

    for f in filenames:
        open = nix.File.open(f, nix.FileMode.ReadOnly) #opens hdf data
        parts = f.split('/')
        data_block = open.blocks[parts[5]] #opens data block of hdf data

        #zugriff auf geschaetzte daten:
        e_pos = data_block.data_arrays['estimated positions']
        estimated_xpos[parts[5]] = e_pos[:,0]  #geschaetzte x positionen
        estimated_ypos[parts[5]] = e_pos[:,1]  #geschaetzte y postitonen
        estimated_pos_times[parts[5]] = e_pos.dimensions[0].ticks #zeitpunkte, an denen geschaetzt wurde
        e_orient = data_block.data_arrays['estimated_orientations']
        estimated_orientations[parts[5]] = e_orient[:,0] #geschaetzte orientierungen - nach oben orientiert = 0 grad

        #zugriff auf tatsaechlich gemessene daten:


        pos = data_block.data_arrays['positions']
        xpos[parts[5]] = (pos[:,0]) #x positionen
        ypos[parts[5]] = (pos[:,1]) #y postitonen
        pos_times[parts[5]] = pos.dimensions[0].ticks
        orient = data_block.data_arrays['orientations']
        orientations[parts[5]] = orient[:,0]

        keys.append(parts[5])


    return estimated_xpos, estimated_ypos, estimated_pos_times, estimated_orientations, xpos, ypos, pos_times, orientations, keys


if __name__ == '__main__':

    E1_coordinates = [485, 60]
    E2_coordinates = [485, 495]
    # 1 Pixel = 0.12 cm

    # funktion, die fuer jeden fisch die passenden h5 filenamen generiert
    h5_filenames1 = get_h5_filenames(videofiles1) #variablen videofiles1-6 stammen aus dem python file 'read_data_versuch4'. videofiles1 enthaelt alle videofilenamen von chip (2015albi02), welche brauchbar sind, videofiles 2, die von chap usw.
    h5_filenames2 = get_h5_filenames(videofiles2) #chap (2015albi01)
    h5_filenames3 = get_h5_filenames(videofiles3) #alfons (2014albi08)
    h5_filenames4 = get_h5_filenames(videofiles4) #trixi (2013albi14)
    h5_filenames5 = get_h5_filenames(videofiles5) #krummschwanz (2013albi09)
    h5_filenames6 = get_h5_filenames(videofiles6) #hermes (2012albi01)


    #funktion soll fuer jeden fisch aus den hdf tracking files die wichtigen variablen wie position, zeit usw auslesen, diese werden dann als dictionaries zurueckgegeben, wobei die filenames als keys dienen
    estimated_xpos1, estimated_ypos1, estimated_pos_times1, estimated_orientations1, xpos1, ypos1, pos_times1, orientations1, keys1 = read_data(h5_filenames1)
    estimated_xpos2, estimated_ypos2, estimated_pos_times2, estimated_orientations2, xpos2, ypos2, pos_times2, orientations2, keys2 = read_data(h5_filenames2)
    estimated_xpos3, estimated_ypos3, estimated_pos_times3, estimated_orientations3, xpos3, ypos3, pos_times3, orientations3, keys3 = read_data(h5_filenames3)
    estimated_xpos4, estimated_ypos4, estimated_pos_times4, estimated_orientations4, xpos4, ypos4, pos_times4, orientations4, keys4 = read_data(h5_filenames4)
    estimated_xpos5, estimated_ypos5, estimated_pos_times5, estimated_orientations5, xpos5, ypos5, pos_times5, orientations5, keys5 = read_data(h5_filenames5)
    estimated_xpos6, estimated_ypos6, estimated_pos_times6, estimated_orientations6, xpos6, ypos6, pos_times6, orientations6, keys6 = read_data(h5_filenames6)

    analyse_tracking_data(xpos1, ypos1, keys1, pos_times1, orientations1, E1_coordinates, E2_coordinates)
    #analyse_tracking_data(xpos2, ypos2, keys2, pos_times2, orientations2, E1_coordinates, E2_coordinates)
    #analyse_tracking_data(xpos3, ypos3, keys3, pos_times3, orientations3, E1_coordinates, E2_coordinates)
    #analyse_tracking_data(xpos4, ypos4, keys4, pos_times4, orientations4, E1_coordinates, E2_coordinates)
    #analyse_tracking_data(xpos5, ypos5, keys5, pos_times5, orientations5, E1_coordinates, E2_coordinates)
    #analyse_tracking_data(xpos6, ypos6, keys6, pos_times6, orientations6, E1_coordinates, E2_coordinates)




