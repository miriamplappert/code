__author__ = 'plappert'
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from IPython import embed
import glob
import nix
from read_data_versuch4 import videofiles1, videofiles2, videofiles3, videofiles4, videofiles5, videofiles6, chosen_electrode1, chosen_electrode2, right_choice_video1, wrong_choice_video1, no_choice_video1, rewarded_electrode_video1, right_choice_video2, wrong_choice_video2, no_choice_video2, rewarded_electrode_video2, chosen_electrode_video1, chosen_electrode_video2
from scipy.stats import linregress
import itertools
from pylab import *
import math
import scipy.io as scio
from compiler.ast import flatten
from matplotlib import gridspec
import scipy.stats as stats


def roc_curve2(E1_positives, E1_false_positives, fish, name):

    E1_true_p = []
    E1_false_p = []


    for i in np.arange(len(E1_positives[0])):
        E1_true_p.append(np.trapz(E1_positives[0][0 : i],dx=0.5))
        E1_false_p.append(np.trapz(E1_false_positives[0][0 : i],dx=0.5))


    roc_auc_E1 = np.trapz(E1_true_p, E1_false_p) #roc auc steht fuer "roc area under curve", es wird also die flaeche unter der roc kurve berechnet
    roc_auc_E1 = "{:2.4}".format(str(roc_auc_E1))
    plt.scatter(E1_false_p, E1_true_p)
    plt.ylabel('True Positiv')
    plt.xlabel('False Positives')
    plt.ylim(0, 1)
    plt.xlim(0, 1)
    plt.title('Elektrode1 ' + name + ' ' + fish)
    plt.text(0.02, 0.95, 'FuK = ' + roc_auc_E1) #FuK steht fuer flaeche unter der kurve
    plt.savefig('roc_curve_E1' + name + fish + '.pdf')
    plt.show()

    return

def distance_histogramm3(E1_distances, E2_distances, rewarded_electrode_video, chosen_electrode_video, fish):

    distances_electrode_was_right_and_chosen = []
    distances_electrode_was_wrong_and_not_chosen = []



    for i in np.arange(len(rewarded_electrode_video)):
        if rewarded_electrode_video[i] == 1 and chosen_electrode_video[i] == 1:
            distances_electrode_was_right_and_chosen.extend(list(E1_distances[i][0, :]))
        elif rewarded_electrode_video[i] == 2 and chosen_electrode_video[i] == 2:
            distances_electrode_was_right_and_chosen.extend(list(E2_distances[i][0, :]))
        elif rewarded_electrode_video[i] == 2 and chosen_electrode_video[i] == 1:
            distances_electrode_was_wrong_and_not_chosen.extend(list(E2_distances[i][0, :]))
        elif rewarded_electrode_video[i] == 1 and chosen_electrode_video[i] == 2:
            distances_electrode_was_wrong_and_not_chosen.extend(list(E1_distances[i][0, :]))


    def use_only_distances_smaler_35cm(distances):
        new_distances = []
        for d in distances:
            if d < 35:
                new_distances.append(d)
        return new_distances

    distances_electrode_was_right_and_chosen = use_only_distances_smaler_35cm(distances_electrode_was_right_and_chosen)
    distances_electrode_was_wrong_and_not_chosen = use_only_distances_smaler_35cm(distances_electrode_was_wrong_and_not_chosen)

    # draw the normalized histogramms:

    f, ((ax1, ax2)) = plt.subplots(2, 1, sharex='col', sharey='row')

    hist1 = ax1.hist(distances_electrode_was_right_and_chosen, np.arange(0, 35.5, 0.5), normed=True)
    hist2 = ax2.hist(distances_electrode_was_wrong_and_not_chosen, np.arange(0, 35.5, 0.5), normed=True)



    ##### Macht den blau karierten Hintergrund der Plots #####
    ax1.grid(color='powderblue', linestyle='-')
    ax1.set_axisbelow(True)
    ax1.spines['bottom'].set_color('powderblue')
    ax1.spines['top'].set_color('powderblue')
    ax1.spines['left'].set_color('powderblue')
    ax1.spines['right'].set_color('powderblue')
    for ticks in ax1.xaxis.get_ticklines() + ax1.yaxis.get_ticklines():
        ticks.set_color('powderblue')

    ax2.grid(color='powderblue', linestyle='-')
    ax2.set_axisbelow(True)
    ax2.spines['bottom'].set_color('powderblue')
    ax2.spines['top'].set_color('powderblue')
    ax2.spines['left'].set_color('powderblue')
    ax2.spines['right'].set_color('powderblue')
    for ticks in ax2.xaxis.get_ticklines() + ax2.yaxis.get_ticklines():
        ticks.set_color('powderblue')

    ###################################################################

    f.canvas.draw()
    plt.savefig('Histogramm_Elektrodendistanzen3_' + fish + '.pdf')
    plt.show()

    return hist1, hist2



def roc_curve(E1_positives, E1_false_positives, E2_positives, E2_false_positives, fish, name):

    E1_true_p = []
    E1_false_p = []
    E2_true_p = []
    E2_false_p = []

    for i in np.arange(len(E1_positives[0])):
        E1_true_p.append(np.trapz(E1_positives[0][0 : i],dx=0.5))
        E1_false_p.append(np.trapz(E1_false_positives[0][0 : i],dx=0.5))
        E2_true_p.append(np.trapz(E2_positives[0][0 : i],dx=0.5))
        E2_false_p.append(np.trapz(E2_false_positives[0][0 : i],dx=0.5))

    roc_auc_E1 = np.trapz(E1_true_p, E1_false_p) #roc auc steht fuer "roc area under curve", es wird also die flaeche unter der roc kurve berechnet
    roc_auc_E1 = "{:2.4}".format(str(roc_auc_E1))
    plt.scatter(E1_false_p, E1_true_p)
    plt.ylabel('True Positiv')
    plt.xlabel('False Positives')
    plt.ylim(0, 1)
    plt.xlim(0, 1)
    plt.title('Elektrode1 ' + name + ' ' + fish)
    plt.text(0.02, 0.95, 'FuK = ' + roc_auc_E1) #FuK steht fuer flaeche unter der kurve
    plt.savefig('roc_curve_E1' + name + fish + '.pdf')
    plt.show()



    roc_auc_E2 = np.trapz(E2_true_p, E2_false_p) #roc auc steht fuer "roc area under curve", es wird also die flaeche unter der roc kurve berechnet
    roc_auc_E2 = "{:2.4}".format(str(roc_auc_E2))
    plt.scatter(E2_false_p, E2_true_p)
    plt.ylabel('True Positiv')
    plt.xlabel('False Positives')
    plt.title('Elektrode2 ' + name + ' ' + fish)
    plt.ylim(0, 1)
    plt.xlim(0, 1)
    plt.text(0.02, 0.95, 'FuK = ' + roc_auc_E2) #FuK steht fuer flaeche unter der kurve
    plt.savefig('roc_curve_E2' + name + fish + '.pdf')
    plt.show()

    return


def velocity_rewarded_and_chosen_electrode(velocities_near_electrodes, velocities_far_electrodes, velocities_near_E1, velocities_near_E2, rewarded_electrode_video, chosen_electrode_video, fish):

    velocity_E1_right_and_chosen = []
    velocity_E1_wrong_and_chosen = []
    velocity_E2_right_and_chosen = []
    velocity_E2_wrong_and_chosen = []

    for i in np.arange(len(rewarded_electrode_video)):
        if rewarded_electrode_video[i] == 1 and chosen_electrode_video[i] == 1:
            if len(velocities_near_E1[i]) > 0:
                velocity_E1_right_and_chosen.extend(list(velocities_near_E1[i][0, :]))
        elif rewarded_electrode_video[i] == 2 and chosen_electrode_video[i] == 1:
             if len(velocities_near_E1[i]) > 0:
                velocity_E1_wrong_and_chosen.extend(list(velocities_near_E1[i][0, :]))
        elif rewarded_electrode_video[i] == 2 and chosen_electrode_video[i] == 2:
             if len(velocities_near_E2[i]) > 0:
                velocity_E2_right_and_chosen.extend(list(velocities_near_E2[i][0, :]))
        elif rewarded_electrode_video[i] == 1 and chosen_electrode_video[i] == 2:
             if len(velocities_near_E2[i]) > 0:
                velocity_E2_wrong_and_chosen.extend(list(velocities_near_E2[i][0, :]))


    f, ((ax, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')

    boxplot_dict = ax.boxplot(velocity_E1_right_and_chosen)
    boxplot_dict2 = ax2.boxplot(velocity_E1_wrong_and_chosen)
    boxplot_dict3 = ax3.boxplot(velocity_E2_right_and_chosen)
    boxplot_dict4 = ax4.boxplot(velocity_E2_wrong_and_chosen)
    ax.set_ylim(0,50)
    ax2.set_ylim(0,50)
    ax3.set_ylim(0,50)
    ax4.set_ylim(0,50)

    ax.set_title('E1 richtig und gewaehlt')
    ax2.set_title('E1 richtig nicht gewaehlt')
    ax3.set_title('E2 richtig und gewaehlt')
    ax4.set_title('E2 richtig nicht gewaehlt')
    ax.set_ylabel("Geschwindigkeit [cm/s]")
    ax2.set_ylabel("Geschwindigkeit [cm/s]")
    ax3.set_ylabel("Geschwindigkeit [cm/s]")
    ax4.set_ylabel("Geschwindigkeit [cm/s]")


    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax4.get_xticklabels(), visible=False)
    plt.savefig('Velocityboxplot2' + fish + '.pdf')
    plt.close()

    print np.median(velocity_E1_right_and_chosen), np.median(velocity_E1_wrong_and_chosen), np.median(velocity_E2_right_and_chosen), np.median(velocity_E2_wrong_and_chosen)

    #u, p_value = stats.mannwhitneyu(velocities_near_electrodes, velocities_far_electrodes, use_continuity=True)
    #print u
    #print p_value

    g, ((x1, x2), (x3, x4)) = plt.subplots(2, 2, sharex='col', sharey='row')

    hist1 = x1.hist(velocity_E1_right_and_chosen, np.arange(0, 50.5, 0.5), normed=True)
    hist2 = x2.hist(velocity_E1_wrong_and_chosen, np.arange(0, 50.5, 0.5), normed=True)
    hist3 = x3.hist(velocity_E2_right_and_chosen, np.arange(0, 50.5, 0.5), normed=True)
    hist4 = x4.hist(velocity_E2_wrong_and_chosen, np.arange(0, 50.5, 0.5), normed=True)

    x1.set_title('E1 richtig und ausgewaehlt')
    x2.set_title('E1 falsch und ausgewaehlt')
    x3.set_title('E2 richtig und ausgewaehlt')
    x4.set_title('E2 falsch und ausgewaehlt')

    ##### Macht den blau karierten Hintergrund der Plots #####
    ax.grid(color='powderblue', linestyle='-')
    ax.set_axisbelow(True)
    ax.spines['bottom'].set_color('powderblue')
    ax.spines['top'].set_color('powderblue')
    ax.spines['left'].set_color('powderblue')
    ax.spines['right'].set_color('powderblue')
    for ticks in ax.xaxis.get_ticklines() + ax.yaxis.get_ticklines():
        ticks.set_color('powderblue')

    ax2.grid(color='powderblue', linestyle='-')
    ax2.set_axisbelow(True)
    ax2.spines['bottom'].set_color('powderblue')
    ax2.spines['top'].set_color('powderblue')
    ax2.spines['left'].set_color('powderblue')
    ax2.spines['right'].set_color('powderblue')
    for ticks in ax2.xaxis.get_ticklines() + ax2.yaxis.get_ticklines():
        ticks.set_color('powderblue')

    ax3.grid(color='powderblue', linestyle='-')
    ax3.set_axisbelow(True)
    ax3.spines['bottom'].set_color('powderblue')
    ax3.spines['top'].set_color('powderblue')
    ax3.spines['left'].set_color('powderblue')
    ax3.spines['right'].set_color('powderblue')
    for ticks in ax3.xaxis.get_ticklines() + ax3.yaxis.get_ticklines():
        ticks.set_color('powderblue')

    ax4.grid(color='powderblue', linestyle='-')
    ax4.set_axisbelow(True)
    ax4.spines['bottom'].set_color('powderblue')
    ax4.spines['top'].set_color('powderblue')
    ax4.spines['left'].set_color('powderblue')
    ax4.spines['right'].set_color('powderblue')
    for ticks in ax4.xaxis.get_ticklines() + ax4.yaxis.get_ticklines():
        ticks.set_color('powderblue')

    ###################################################################

    g.canvas.draw()
    plt.savefig('Histogramm_Geschwindigkeiten_mit_Fischentscheidung' + fish + '.pdf')
    plt.close()

    return hist1, hist2, hist3, hist4


def distances_electrodes_histogramm_with_fish_choice(E1_distances, E2_distances, rewarded_electrode_video, chosen_electrode_video, fish):


    E1_distances_1_was_right_and_chosen = []
    E1_distances_1_was_wrong_and_chosen = []
    E2_distances_2_was_right_and_chosen = []
    E2_distances_2_was_wrong_and_chosen = []


    for i in np.arange(len(rewarded_electrode_video)):
        if rewarded_electrode_video[i] == 1 and chosen_electrode_video[i] == 1:
            E1_distances_1_was_right_and_chosen.extend(list(E1_distances[i][0, :]))
        elif rewarded_electrode_video[i] == 2 and chosen_electrode_video[i] == 1:
            E1_distances_1_was_wrong_and_chosen.extend(list(E1_distances[i][0, :]))
        elif rewarded_electrode_video[i] == 2 and chosen_electrode_video[i] == 2:
            E2_distances_2_was_right_and_chosen.extend(list(E2_distances[i][0, :]))
        elif rewarded_electrode_video[i] == 1 and chosen_electrode_video[i] == 2:
            E2_distances_2_was_wrong_and_chosen.extend(list(E2_distances[i][0, :]))



    def use_only_distances_smaler_35cm(distances):
        new_distances = []
        for d in distances:
            if d < 35:
                new_distances.append(d)
        return new_distances

    E1_distances_1_was_right_and_chosen = use_only_distances_smaler_35cm(E1_distances_1_was_right_and_chosen)
    E1_distances_1_was_wrong_and_chosen = use_only_distances_smaler_35cm(E1_distances_1_was_wrong_and_chosen)
    E2_distances_2_was_right_and_chosen = use_only_distances_smaler_35cm(E2_distances_2_was_right_and_chosen)
    E2_distances_2_was_wrong_and_chosen = use_only_distances_smaler_35cm(E2_distances_2_was_wrong_and_chosen)


    # draw the normalized histogramms:

    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')

    hist1 = ax1.hist(E1_distances_1_was_right_and_chosen, np.arange(0, 35.5, 0.5), normed=True)
    hist2 = ax2.hist(E1_distances_1_was_wrong_and_chosen, np.arange(0, 35.5, 0.5), normed=True)
    hist3 = ax3.hist(E2_distances_2_was_right_and_chosen, np.arange(0, 35.5, 0.5), normed=True)
    hist4 = ax4.hist(E2_distances_2_was_wrong_and_chosen, np.arange(0, 35.5, 0.5), normed=True)

    ax1.set_title('E1 richtig positiv')
    ax2.set_title('E1 falsch positiv')
    ax3.set_title('E2 richtig positiv')
    ax4.set_title('E2 falsch positiv')

    ax1.set_ylabel('Haeufigkeit')
    ax3.set_ylabel('Haeufigkeit')
    ax3.set_xlabel('Distanz [Zentimeter]')
    ax4.set_xlabel('Distanz [Zentimeter]')

    if 'albi02' in fish:
        ax1.set_ylim(0,.14)
        ax2.set_ylim(0,.14)
        ax3.set_ylim(0,.14)
        ax4.set_ylim(0,.14)

    if 'albi01' in fish:
        ax1.set_ylim(0,.25)
        ax2.set_ylim(0,.25)
        ax3.set_ylim(0,.25)
        ax4.set_ylim(0,.25)

    ##### Macht den blau karierten Hintergrund der Plots #####
    ax1.grid(color='powderblue', linestyle='-')
    ax1.set_axisbelow(True)
    ax1.spines['bottom'].set_color('powderblue')
    ax1.spines['top'].set_color('powderblue')
    ax1.spines['left'].set_color('powderblue')
    ax1.spines['right'].set_color('powderblue')
    for ticks in ax1.xaxis.get_ticklines() + ax1.yaxis.get_ticklines():
        ticks.set_color('powderblue')

    ax2.grid(color='powderblue', linestyle='-')
    ax2.set_axisbelow(True)
    ax2.spines['bottom'].set_color('powderblue')
    ax2.spines['top'].set_color('powderblue')
    ax2.spines['left'].set_color('powderblue')
    ax2.spines['right'].set_color('powderblue')
    for ticks in ax2.xaxis.get_ticklines() + ax2.yaxis.get_ticklines():
        ticks.set_color('powderblue')

    ax3.grid(color='powderblue', linestyle='-')
    ax3.set_axisbelow(True)
    ax3.spines['bottom'].set_color('powderblue')
    ax3.spines['top'].set_color('powderblue')
    ax3.spines['left'].set_color('powderblue')
    ax3.spines['right'].set_color('powderblue')
    for ticks in ax3.xaxis.get_ticklines() + ax3.yaxis.get_ticklines():
        ticks.set_color('powderblue')

    ax4.grid(color='powderblue', linestyle='-')
    ax4.set_axisbelow(True)
    ax4.spines['bottom'].set_color('powderblue')
    ax4.spines['top'].set_color('powderblue')
    ax4.spines['left'].set_color('powderblue')
    ax4.spines['right'].set_color('powderblue')
    for ticks in ax4.xaxis.get_ticklines() + ax4.yaxis.get_ticklines():
        ticks.set_color('powderblue')

    ###################################################################

    f.canvas.draw()
    plt.savefig('Histogramm_Elektrodendistanzen_mit_Fischentscheidung' + fish + '.pdf')
    plt.close()

    return hist1, hist2, hist3, hist4


def distances_electrodes_histogramm(E1_distances, E2_distances, rewarded_electrode_video, fish):
    E1_distances_1_was_right = []
    E2_distances_1_was_right = []
    E1_distances_2_was_right = []
    E2_distances_2_was_right = []


    for i in np.arange(len(rewarded_electrode_video)):
        if rewarded_electrode_video[i] == 1:
            E1_distances_1_was_right.extend(list(E1_distances[i][0, :]))
            E2_distances_1_was_right.extend(list(E2_distances[i][0, :]))
        elif rewarded_electrode_video[i] == 2:
            E1_distances_2_was_right.extend(list(E1_distances[i][0, :]))
            E2_distances_2_was_right.extend(list(E2_distances[i][0, :]))



    def use_only_distances_smaler_35cm(distances):
        new_distances = []
        for d in distances:
            if d < 35:
                new_distances.append(d)
        return new_distances

    E1_distances_1_was_right = use_only_distances_smaler_35cm(E1_distances_1_was_right)
    E2_distances_1_was_right = use_only_distances_smaler_35cm(E2_distances_1_was_right)
    E1_distances_2_was_right = use_only_distances_smaler_35cm(E1_distances_2_was_right)
    E2_distances_2_was_right = use_only_distances_smaler_35cm(E2_distances_2_was_right)


    # draw the normalized histogramms:

    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')

    hist1 = ax1.hist(E1_distances_1_was_right, 50, normed=True)
    hist2 = ax2.hist(E1_distances_2_was_right, 50, normed=True)
    hist3 = ax3.hist(E2_distances_2_was_right, 50, normed=True)
    hist4 = ax4.hist(E2_distances_1_was_right, 50, normed=True)


    ax1.set_title('E1 richtige')
    ax2.set_title('E1 falsche')
    ax3.set_title('E2 richtige')
    ax4.set_title('E2 falsche')

    f.canvas.draw()
    plt.savefig('Histogramm_Elektrodendistanzen' + fish + '.pdf')
    plt.close()


    return



def velocity_box_plot(velocities_near_electrodes, velocities_far_electrodes, velocities_near_E1, velocities_near_E2, fish):


    velocities_near_electrodes_biglist = []
    velocities_far_electrodes_biglist = []
    velocities_near_E1_biglist = []
    velocities_near_E2_biglist = []

    for i in np.arange(len(velocities_near_electrodes)):
        if len(velocities_near_electrodes[i]) > 0:
            velocities_near_electrodes_biglist.extend(list(velocities_near_electrodes[i][0, :]))
        if len(velocities_far_electrodes[i]) > 0:
            velocities_far_electrodes_biglist.extend(list(velocities_far_electrodes[i][0, :]))
        if len(velocities_near_E1[i]) > 0:
            velocities_near_E1_biglist.extend(list(velocities_near_E1[i][0, :]))
        if len(velocities_near_E2[i]) > 0:
            velocities_near_E2_biglist.extend(list(velocities_near_E2[i][0, :]))



    f, ((ax, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')

    boxplot_dict = ax.boxplot(velocities_near_electrodes_biglist)
    boxplot_dict2 = ax2.boxplot(velocities_far_electrodes_biglist)
    boxplot_dict3 = ax3.boxplot(velocities_near_E1_biglist)
    boxplot_dict4 = ax4.boxplot(velocities_near_E2_biglist)

    ax.set_ylim(0,50)
    ax2.set_ylim(0,50)
    ax3.set_ylim(0,50)
    ax4.set_ylim(0,50)

    ax.set_title('In der Naehe der Elektroden')
    ax2.set_title('Entfernt von den Elektroden')
    ax3.set_title('Nahe E1')
    ax4.set_title('Nahe E2')
    ax.set_ylabel("Geschwindigkeit [cm/s]")
    ax2.set_ylabel("Geschwindigkeit [cm/s]")
    ax3.set_ylabel("Geschwindigkeit [cm/s]")
    ax4.set_ylabel("Geschwindigkeit [cm/s]")



        ##### Macht den blau karierten Hintergrund der Plots #####
    ax.grid(color='powderblue', linestyle='-')
    ax.set_axisbelow(True)
    ax.spines['bottom'].set_color('powderblue')
    ax.spines['top'].set_color('powderblue')
    ax.spines['left'].set_color('powderblue')
    ax.spines['right'].set_color('powderblue')
    for ticks in ax.xaxis.get_ticklines() + ax.yaxis.get_ticklines():
        ticks.set_color('powderblue')

    ax2.grid(color='powderblue', linestyle='-')
    ax2.set_axisbelow(True)
    ax2.spines['bottom'].set_color('powderblue')
    ax2.spines['top'].set_color('powderblue')
    ax2.spines['left'].set_color('powderblue')
    ax2.spines['right'].set_color('powderblue')
    for ticks in ax2.xaxis.get_ticklines() + ax2.yaxis.get_ticklines():
        ticks.set_color('powderblue')

    ax3.grid(color='powderblue', linestyle='-')
    ax3.set_axisbelow(True)
    ax3.spines['bottom'].set_color('powderblue')
    ax3.spines['top'].set_color('powderblue')
    ax3.spines['left'].set_color('powderblue')
    ax3.spines['right'].set_color('powderblue')
    for ticks in ax3.xaxis.get_ticklines() + ax3.yaxis.get_ticklines():
        ticks.set_color('powderblue')

    ax4.grid(color='powderblue', linestyle='-')
    ax4.set_axisbelow(True)
    ax4.spines['bottom'].set_color('powderblue')
    ax4.spines['top'].set_color('powderblue')
    ax4.spines['left'].set_color('powderblue')
    ax4.spines['right'].set_color('powderblue')
    for ticks in ax4.xaxis.get_ticklines() + ax4.yaxis.get_ticklines():
        ticks.set_color('powderblue')

    ###################################################################

    '''
    for b in boxplot_dict['fliers']:
        b.set_color('white')
    for b in boxplot_dict2['fliers']:
        b.set_color('white')
    '''
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax4.get_xticklabels(), visible=False)
    plt.savefig('Velocityboxplot' + fish + '.pdf')
    plt.close()


    #u, p_value = stats.mannwhitneyu(velocities_near_electrodes, velocities_far_electrodes, use_continuity=True)
    #print u
    #print p_value


def get_mismatch_filenames(mismatch_indices, list_of_filenames, estimated_decision):
    mismatch_filenames = []
    computer_decision = []
    for i in np.arange(len(mismatch_indices)):
        mismatch_filenames.append(list_of_filenames[i])
        computer_decision.append(estimated_decision[i])

    return mismatch_filenames, computer_decision

def compare_estimated_to_real_decision(estimated_decision, chosen_electrode):

    estimated_decision = estimated_decision.tolist()

    print estimated_decision

    mismatch = []

    for i in np.arange(len(chosen_electrode)): # for loop goes through the indices of the chosen_electrodes list
        if estimated_decision[i] == chosen_electrode[i]:
            continue
        else:
            mismatch.append(i)


    return mismatch

def decision_maker(small_E1_distance_amount, small_E2_distance_amount, orientation_divergence_E1, orientation_divergence_E2, k):

    if small_E2_distance_amount == 0.0 and small_E1_distance_amount == 0.0:
        estimated_decision_after_distance = 0
    elif small_E2_distance_amount == 0.0:
        estimated_decision_after_distance = 1
    elif small_E1_distance_amount == 0.0:
        estimated_decision_after_distance = 2
    else:
        E1_E2_distance_proportion = float(small_E1_distance_amount) / float(small_E2_distance_amount)

        a = E1_E2_distance_proportion

        if a <= 1/3:
            estimated_decision_after_distance = 2
        if a >= 3:
            estimated_decision_after_distance = 1
        else:
            estimated_decision_after_distance = 0

    '''
    if np.mean(orientation_divergence_E1) < np.mean(orientation_divergence_E2):
        estimated_decision_after_orientation = 1
    if np.mean(orientation_divergence_E2) < np.mean(orientation_divergence_E1):
        estimated_decision_after_orientation = 2
    else:
        estimated_decision_after_orientation = 0


    print estimated_decision_after_distance, estimated_decision_after_orientation
    '''

    return estimated_decision_after_distance


def velocity_near_and_far_from_electrodes(velocities,distance_to_E1, distance_to_E2, filename):

    velocities_near_electrodes = []
    velocities_far_from_electrodes = []
    velocities_near_E1 = []
    velocities_near_E2 = []

    for c in np.arange(len(distance_to_E1)): # for schleife erstellt liste mit den indices, an denen im eod array eine 0 steht

        if velocities[c-1] > 150:
            continue
        elif distance_to_E1[c] < (100*0.12) or distance_to_E2 < (100*0.12): # 100*12 fuer die umrechnung von pixeln in zentimeter --> 100 Pixel werden durch den faktor 0,12 in zentimeter umgerechnet
            velocities_near_electrodes.append(velocities[c-1])
        else:
            velocities_far_from_electrodes.append(velocities[c-1])


    for d in np.arange(len(distance_to_E1)):
        if distance_to_E1[d] < (100*0.12):
            velocities_near_E1.append(velocities[d-1])
        if distance_to_E2[d] < (100*0.12):
            velocities_near_E2.append(velocities[d-1])


    return velocities_near_electrodes, velocities_far_from_electrodes, velocities_near_E1, velocities_near_E2

def orientation_near_electrode_plot(x_position_near_e1, y_position_near_e1, orientation_near_e1, x_position_near_e2, y_position_near_e2, orientation_near_e2, E1_coordinates, E2_coordinates, velocities_E1, velocities_E2, x_pos, y_pos, k):
    print '1', len (velocities_E2)
    print '1', len(orientation_near_e2)


    #only velocities < 100:
    velocity_E1 = []
    velocity_E2 = []

    if len(velocities_E1) > 0:
        for v in np. arange(len(velocities_E1)):
            if velocities_E1[v] > 100:
                velocity_E1.append(0)
            else:
                velocity_E1.append(velocities_E1[v])

    if len(velocities_E2) > 0:
        for v2 in np. arange(len(velocities_E2)):
            if velocities_E2[v2] > 100:
                velocity_E2.append(0)
            else:
                velocity_E2.append(velocities_E2[v2])

    #normalize velocities:
    normalized_E1_velocities = []
    normalized_E2_velocities = []

    if len(velocity_E1) > 0:
        for i1 in np.arange(len(velocity_E1)):
            v_neu = (velocity_E1[i1] - min(velocity_E1)) / (max(velocity_E1) - min(velocity_E1))
            normalized_E1_velocities.append(v_neu)

    if len(velocity_E2) > 0:
        for i2 in np.arange(len(velocity_E2)):
            v_neu2 = (velocity_E2[i2] - min(velocity_E2)) / (max(velocity_E2) - min(velocity_E2))
            normalized_E2_velocities.append(v_neu2)

    img = plt.imread('/home/plappert/Data/video_files/' + k + '/' + k + '_OV_path.png')

    u = []
    v = []
    x = x_position_near_e1
    y = y_position_near_e1

    for i in np.arange(len(orientation_near_e1)):
        theta = (360 - float(orientation_near_e1[i])) - 90
        theta = pi/180 * theta # umrechnung von grad in bogenmass
        if len(normalized_E1_velocities) > 0:
            r = 5 * normalized_E1_velocities[i]; # magnitude (length) of arrow
        else:
            r = 1
        u.append(r * np.cos(theta)) # convert polar (theta,r) to cartesian
        v.append(r * np.sin(theta))

    h = plt.quiver(x,y,u,v, color = 'green', scale_units='xy', scale=0.1)
    plt.grid()
    kreis_electrode1 = plt.Circle(E1_coordinates,5,color='r')

    print len(normalized_E2_velocities)
    print len(orientation_near_e2)

    u2 = []
    v2 = []
    x2 = x_position_near_e2
    y2 = y_position_near_e2



    for j in np.arange(len(orientation_near_e2)):
         beta = (360 - float(orientation_near_e2[j])) - 90 # umrechnen der vom trackingprogramm gespeicherten orientierung, bei der eine suedliche ausrichtung null grad entspricht und anschliesend im uhrzeigersinn gedreht wird, in die ausrichtung die von quiver benutzt wird, bei der null grad oestlicher ausrichtung entspricht und gegen den uhrzeigersinn gedreht wird
         beta = pi/180 * beta # umrechnung von grad in bogenmass
         if len(normalized_E2_velocities) > 0:
            r = 5 * normalized_E2_velocities[j]  # magnitude (length) of arrow to plot
         else:
            r = 1
         u2.append(r * np.cos(beta)) # convert polar (theta,r) to cartesian
         v2.append(r * np.sin(beta))

    h = plt.quiver(x2,y2,u2,v2, color='green', scale_units='xy', scale=0.1)
    plt.grid()
    plt.imshow(img, origin='lower')
    kreis_electrode2 = plt.Circle(E2_coordinates,5,color='b')
    fig = plt.gcf()
    fig.gca().add_artist(kreis_electrode1)
    fig.gca().add_artist(kreis_electrode2)
    plt.ylim(50, 510)
    plt.xlim(50, 700)
    plt.savefig('Orientierung_' + k + '.pdf')
    plt.show()


def orientation_to_electrode(orientations, x_positions, y_positions, distance_to_E1, distance_to_E2, E1_coordinates, E2_coordinates, k):

    index_great_distance_E1 = []
    for i in np.arange(len(distance_to_E1)): # for schleife erstellt liste mit den indices, an denen im eod array eine 0 steht
        if distance_to_E1[i] > (100*0.12): # 100*12 fuer die umrechnung von pixeln in zentimeter --> 100 Pixel werden durch den faktor 0,12 in zentimeter umgerechnet
            index_great_distance_E1.append(i)

    index_great_distance_E2 = []
    for i2 in np.arange(len(distance_to_E2)): # for schleife erstellt liste mit den indices, an denen im eod array eine 0 steht
        if distance_to_E2[i2] > (100*0.12):
            index_great_distance_E2.append(i2)


    x_position_near_e1 = [i for j, i in enumerate(x_positions) if j not in index_great_distance_E1]
    y_position_near_e1 = [i for j, i in enumerate(y_positions) if j not in index_great_distance_E1]
    orientation_near_e1 = [i for j, i in enumerate(orientations) if j not in index_great_distance_E1]

    x_position_near_e2 = [i2 for j2, i2 in enumerate(x_positions) if j2 not in index_great_distance_E2]
    y_position_near_e2 = [i2 for j2, i2 in enumerate(y_positions) if j2 not in index_great_distance_E2]
    orientation_near_e2 = [i2 for j2, i2 in enumerate(orientations) if j2 not in index_great_distance_E2]


    theoretical_orientation_E1 = []

    for f in np.arange(len(x_position_near_e1)):
        gegenkathete = np.abs(float(x_position_near_e1[f]) - E1_coordinates[0])
        hypothenuse = np.sqrt(((x_position_near_e1[f] - E1_coordinates[0])**2) + ((y_position_near_e1[f] - E1_coordinates[1])**2))
        a = gegenkathete / hypothenuse #sin(alpha) = gegenkathete/hypothenuse = Betrag von (xpositionfisch - xpositionelektrode1) / distanz zu elektrode 1
        alpha = math.asin(a)
        if x_position_near_e1[f] < E1_coordinates[0]:
            alpha = 360 - alpha
            theoretical_orientation_E1.append(alpha)
        if x_position_near_e1[f] >= E1_coordinates[0]:
            alpha = alpha
            theoretical_orientation_E1.append(alpha)


    theoretical_orientation_E2 = []

    for e in np.arange(len(x_position_near_e2)):
        gegenkathete2 = np.abs(float(x_position_near_e2[e]) - E2_coordinates[0])
        hypothenuse2 = np.sqrt(((x_position_near_e2[e] - E2_coordinates[0])**2) + ((y_position_near_e2[e] - E2_coordinates[1])**2))
        a2 = gegenkathete2 / hypothenuse2 #sin(alpha) = gegenkathete/hypothenuse = Betrag von (xpositionfisch - xpositionelektrode1) / distanz zu elektrode 1
        alpha2 = math.asin(a2)
        if x_position_near_e2[e] < E2_coordinates[0]:
            alpha2 = 180 - alpha2
            theoretical_orientation_E2.append(alpha2)
        if x_position_near_e2[e] >= E2_coordinates[0]:
            alpha2 = 180 + alpha2
            theoretical_orientation_E2.append(alpha2)

    orientation_divergence_E1 = []
    orientation_divergence_E2 = []

    for n in np.arange(len(theoretical_orientation_E1)):
        orientation_divergence_E1.append(np.abs(theoretical_orientation_E1[n] - orientation_near_e1[n]))
    for m in np.arange(len(theoretical_orientation_E2)):
        orientation_divergence_E2.append(np.abs(theoretical_orientation_E2[m] - orientation_near_e2[m]))


    return x_position_near_e1, y_position_near_e1, orientation_near_e1, x_position_near_e2, y_position_near_e2, orientation_near_e2, orientation_divergence_E1, orientation_divergence_E2


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
    estimated_decision = []
    big_ls_vels_near_elects = []
    big_ls_vels_far_elects = []
    big_ls_vels_near_e1 = []
    big_ls_vels_near_e2 = []
    E1_distances = []
    E2_distances = []
    big_ls_velocities = []

    for k in keys: #for loop goes through the filenames
        x_positions = xpos[k] # calls list under the respective filename k
        y_positions = ypos[k]
        pos_times = pos_time[k]
        orientations = orientation[k]

        #further analysis functions:

        distance_to_E1, distance_to_E2 = get_distance_to_electrode(x_positions, y_positions, E1_coordinates, E2_coordinates)
        #distance_time_plot(distance_to_E1, distance_to_E2, pos_times, k)
        velocities = get_velocity(x_positions, y_positions, pos_times)
        #velocity_time_plot(velocities, pos_times, k)
        small_E1_distance_amount, small_E2_distance_amount = small_distance_to_electrode_abundance(distance_to_E1, distance_to_E2)
        #small_distance_amount_bar_plot(small_E1_distance_amount, small_E2_distance_amount, k)
        #distance_velocity_plot(distance_to_E1, distance_to_E2, velocities, k)
        x_position_near_e1, y_position_near_e1, orientation_near_e1, x_position_near_e2, y_position_near_e2, orientation_near_e2, orientation_divergence_E1, orientation_divergence_E2 = orientation_to_electrode(orientations, x_positions, y_positions, distance_to_E1, distance_to_E2, E1_coordinates, E2_coordinates, k)
        velocities_near_electrodes, velocities_far_from_electrodes, velocities_near_E1,  velocities_near_E2 = velocity_near_and_far_from_electrodes(velocities, distance_to_E1, distance_to_E2, k)
        orientation_near_electrode_plot(x_position_near_e1, y_position_near_e1, orientation_near_e1, x_position_near_e2, y_position_near_e2, orientation_near_e2, E1_coordinates, E2_coordinates, velocities_near_E1, velocities_near_E2, x_positions, y_positions, k)
        estimated_decision_after_distance = decision_maker(small_E1_distance_amount, small_E2_distance_amount, orientation_divergence_E1, orientation_divergence_E2, k)


        big_ls_vels_near_e1.append(velocities_near_E1)
        big_ls_vels_near_e2.append(velocities_near_E2)
        big_ls_vels_near_elects.append(velocities_near_electrodes)
        big_ls_vels_far_elects.append(velocities_far_from_electrodes)
        estimated_decision.append(estimated_decision_after_distance)
        E1_distances.append(distance_to_E1)
        E2_distances.append(distance_to_E2)
        big_ls_velocities.append(velocities)

    #big_ls_vels_near_elects = flatten(big_ls_vels_near_elects)
    #big_ls_vels_far_elects = flatten(big_ls_vels_far_elects)


    return estimated_decision, big_ls_vels_near_elects, big_ls_vels_far_elects, E1_distances, E2_distances, big_ls_vels_near_e1, big_ls_vels_near_e2, big_ls_velocities
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

    E1_coordinates = [485, 58]
    E2_coordinates = [485, 495]
    # 1 Pixel = 0.12 cm

    fish1 = '2015albi02'
    fish2 = '2015albi01'
    fish3 = '2014albi08'
    fish4 = '2013albi14'
    fish5 = '2013albi09'
    fish6 = '2012albi01'

    saved_data = glob.glob('analysed_data.mat')
    saved_data2 = glob.glob('velocities.mat')
    if len(saved_data) > 100000000000:
        analysed_data = scio.loadmat('analysed_data.mat')
        velocities = scio.loadmat('velocities.mat')

        estimated_decision1 = analysed_data['estimated_decision1']
        estimated_decision2 = analysed_data['estimated_decision2']
        estimated_decision1 = estimated_decision1[0]
        estimated_decision2 = estimated_decision2[0]
        keys1 = analysed_data['keys1']
        keys2 = analysed_data['keys2']

        velocities_near_electrodes1 = velocities['velocities_near_electrodes1']
        velocities_far_electrodes1 = velocities['velocities_far_electrodes1']
        velocities_near_electrodes2 = velocities['velocities_near_electrodes2']
        velocities_far_electrodes2 = velocities['velocities_far_electrodes2']
        velocities_near_electrodes1 = velocities_near_electrodes1[0]
        velocities_far_electrodes1 = velocities_far_electrodes1[0]
        velocities_near_electrodes2 =velocities_near_electrodes2[0]
        velocities_far_electrodes2 = velocities_far_electrodes2[0]

        velocities_near_E1_1 = velocities['velocities_near_E1_1']
        velocities_near_E2_1 = velocities['velocities_near_E2_1']
        velocities_near_E1_2 = velocities['velocities_near_E1_2']
        velocities_near_E2_2 = velocities['velocities_near_E2_2']
        velocities_near_E1_1 = velocities_near_E1_1[0]
        velocities_near_E2_1 = velocities_near_E2_1[0]
        velocities_near_E1_2 = velocities_near_E1_2[0]
        velocities_near_E2_2 = velocities_near_E2_2[0]

        velocities1 = velocities['velocities1']
        velocities2 = velocities['velocities2']
        velocities1 = velocities1[0]
        velocities2 = velocities2[0]

        E1_distances1 = analysed_data['E1_distances1']
        E2_distances1 = analysed_data['E2_distances1']
        E1_distances2 = analysed_data['E1_distances2']
        E2_distances2 = analysed_data['E2_distances2']
        E1_distances1 = E1_distances1[0]
        E2_distances1 = E2_distances1[0]
        E1_distances2 = E1_distances2[0]
        E2_distances2 = E2_distances2[0]

        #distances_electrodes_histogramm(E1_distances1, E2_distances1, rewarded_electrode_video1, fish1)
        #distances_electrodes_histogramm(E1_distances2, E2_distances2, rewarded_electrode_video2, fish2)

        E1_postives_distance1, E1_false_positives_distance1, E2_positives_distance1, E2_false_positives_distance1 = distances_electrodes_histogramm_with_fish_choice(E1_distances1, E2_distances1, rewarded_electrode_video1, chosen_electrode_video1, fish1)
        E1_postives_distance2, E1_false_positives_distance2, E2_positives_distance2, E2_false_positives_distance2 = distances_electrodes_histogramm_with_fish_choice(E1_distances2, E2_distances2, rewarded_electrode_video2, chosen_electrode_video2, fish2)

        #roc_curve(E1_postives_distance1, E1_false_positives_distance1, E2_positives_distance1, E2_false_positives_distance1, fish1, 'Distanz')
        #roc_curve(E1_postives_distance2, E1_false_positives_distance2, E2_positives_distance2, E2_false_positives_distance2, fish2, 'Distanz')

        hist_distances_right_and_chosen1, hist_distances_wrong_and_not_chosen1 = distance_histogramm3(E1_distances1, E2_distances1, rewarded_electrode_video1, chosen_electrode_video1, fish1)
        hist_distances_right_and_chosen2, hist_distances_wrong_and_not_chosen2 = distance_histogramm3(E1_distances2, E2_distances2, rewarded_electrode_video2, chosen_electrode_video2, fish2)

        roc_curve2(hist_distances_right_and_chosen1, hist_distances_wrong_and_not_chosen1, fish1, 'Distanz2')
        roc_curve2(hist_distances_right_and_chosen2, hist_distances_wrong_and_not_chosen2, fish2, 'Distanz2')

        '''
        mismatch_indices1 = compare_estimated_to_real_decision(estimated_decision1, chosen_electrode1)
        mismatch_indices2 = compare_estimated_to_real_decision(estimated_decision2, chosen_electrode2)

        mismatch_filenames1, computer_decision1 = get_mismatch_filenames(mismatch_indices1, keys1, estimated_decision1)
        mismatch_filenames2, computer_decision2 = get_mismatch_filenames(mismatch_indices2, keys2, estimated_decision2)

        print mismatch_filenames1, computer_decision1
        print mismatch_filenames2, computer_decision2

        '''

        velocity_box_plot(velocities_near_electrodes1, velocities_far_electrodes1, velocities_near_E1_1, velocities_near_E2_1, fish1)
        velocity_box_plot(velocities_near_electrodes2, velocities_far_electrodes2,velocities_near_E1_2, velocities_near_E2_2, fish2)

        E1_postives_velocity1, E1_false_positives_velocity1, E2_positives_velocity1, E2_false_positives_velocity1 = velocity_rewarded_and_chosen_electrode(velocities_near_electrodes1, velocities_far_electrodes1, velocities_near_E1_1, velocities_near_E2_1, rewarded_electrode_video1, chosen_electrode_video1,fish1)
        E1_postives_velocity2, E1_false_positives_velocity2, E2_positives_velocity2, E2_false_positives_velocity2 = velocity_rewarded_and_chosen_electrode(velocities_near_electrodes2, velocities_far_electrodes2, velocities_near_E1_2, velocities_near_E2_2, rewarded_electrode_video2, chosen_electrode_video2,fish2)

        #roc_curve(E1_postives_velocity1, E1_false_positives_velocity1, E2_positives_velocity1, E2_false_positives_velocity1, fish1, 'Geschwindigkeit')
        #roc_curve(E1_postives_velocity2, E1_false_positives_velocity2, E2_positives_velocity2, E2_false_positives_velocity2, fish2, 'Geschwindigkeit')


    else:

        # funktion, die fuer jeden fisch die passenden h5 filenamen generiert
        h5_filenames1 = get_h5_filenames(videofiles1) #variablen videofiles1-6 stammen aus dem python file 'read_data_versuch4'. videofiles1 enthaelt alle videofilenamen von chip (2015albi02), welche brauchbar sind, videofiles 2, die von chap usw.
        h5_filenames2 = get_h5_filenames(videofiles2) #chap (2015albi01)
        #h5_filenames3 = get_h5_filenames(videofiles3) #alfons (2014albi08)
        #h5_filenames4 = get_h5_filenames(videofiles4) #trixi (2013albi14)
        #h5_filenames5 = get_h5_filenames(videofiles5) #krummschwanz (2013albi09)
        #h5_filenames6 = get_h5_filenames(videofiles6) #hermes (2012albi01)


        #funktion soll fuer jeden fisch aus den hdf tracking files die wichtigen variablen wie position, zeit usw auslesen, diese werden dann als dictionaries zurueckgegeben, wobei die filenames als keys dienen
        estimated_xpos1, estimated_ypos1, estimated_pos_times1, estimated_orientations1, xpos1, ypos1, pos_times1, orientations1, keys1 = read_data(h5_filenames1)
        estimated_xpos2, estimated_ypos2, estimated_pos_times2, estimated_orientations2, xpos2, ypos2, pos_times2, orientations2, keys2 = read_data(h5_filenames2)
        #estimated_xpos3, estimated_ypos3, estimated_pos_times3, estimated_orientations3, xpos3, ypos3, pos_times3, orientations3, keys3 = read_data(h5_filenames3)
        #estimated_xpos4, estimated_ypos4, estimated_pos_times4, estimated_orientations4, xpos4, ypos4, pos_times4, orientations4, keys4 = read_data(h5_filenames4)
        #estimated_xpos5, estimated_ypos5, estimated_pos_times5, estimated_orientations5, xpos5, ypos5, pos_times5, orientations5, keys5 = read_data(h5_filenames5)
        #estimated_xpos6, estimated_ypos6, estimated_pos_times6, estimated_orientations6, xpos6, ypos6, pos_times6, orientations6, keys6 = read_data(h5_filenames6)

        estimated_decision1, velocities_near_electrodes1, velocities_far_electrodes1, E1_distances1, E2_distances1, velocities_near_E1_1, velocities_nearE2_1, velocities1 = analyse_tracking_data(xpos1, ypos1, keys1, pos_times1, orientations1, E1_coordinates, E2_coordinates)
        estimated_decision2, velocities_near_electrodes2, velocities_far_electrodes2, E1_distances2, E2_distances2, velocities_near_E1_2, velocities_nearE2_2, velocities2 = analyse_tracking_data(xpos2, ypos2, keys2, pos_times2, orientations2, E1_coordinates, E2_coordinates)
        #estimated_decision3, velocities_near_electrodes3, velocities_far_electrodes3 = analyse_tracking_data(xpos3, ypos3, keys3, pos_times3, orientations3, E1_coordinates, E2_coordinates)
        #estimated_decision4, velocities_near_electrodes4, velocities_far_electrodes4 = analyse_tracking_data(xpos4, ypos4, keys4, pos_times4, orientations4, E1_coordinates, E2_coordinates)
        #estimated_decision5, velocities_near_electrodes5, velocities_far_electrodes5 = analyse_tracking_data(xpos5, ypos5, keys5, pos_times5, orientations5, E1_coordinates, E2_coordinates)
        #estimated_decision6, velocities_near_electrodes6, velocities_far_electrodes6 = analyse_tracking_data(xpos6, ypos6, keys6, pos_times6, orientations6, E1_coordinates, E2_coordinates)

        scio.savemat('analysed_data.mat', {'estimated_decision1': estimated_decision1,
                                           'estimated_decision2':estimated_decision2,
                                           'keys1':keys1, 'keys2':keys2,
                                           'E1_distances1' : E1_distances1,
                                           'E2_distances1': E2_distances1,
                                           'E1_distances2' : E1_distances2,
                                           'E2_distances2' : E2_distances2})


        scio.savemat('velocities.mat', {'velocities_near_electrodes1': velocities_near_electrodes1,
                                        'velocities_far_electrodes1': velocities_far_electrodes1,
                                        'velocities_near_electrodes2': velocities_near_electrodes2,
                                        'velocities_far_electrodes2': velocities_far_electrodes2,
                                        'velocities_near_E1_1':velocities_near_E1_1,
                                        'velocities_near_E2_1': velocities_nearE2_1,
                                        'velocities_near_E1_2':velocities_near_E1_2,
                                        'velocities_near_E2_2': velocities_nearE2_2,
                                        'velocities1': velocities1,
                                        'velocities2':velocities2})
