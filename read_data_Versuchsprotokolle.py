#! /usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'plappert'
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from IPython import embed
import glob
from scipy.stats import linregress
import plotly.plotly as py
from pylab import *
import seaborn as sns
from compiler.ast import flatten

def date_eod_temperature_plot(dates, eods, temperatures, fish):
    """
    function that creats a plot with two y-axis to show the correlation between temperature and eod over the time

    :param dates: array of all dates from all trials of the fish
    :param eods: array of all eods from all trials of the fish
    :param temperatures: array of all temperatures from all trials of the fish
    :param fish: contains fish id (e.g.: 2015albi02)
    :param enu: contains number of file
    :return: plot with two y-axis: first y-axis=mean_eod, second y_axis=mean_temperature, x-axis=date of trial
    """


    index_is_zero = []
    for i in np.arange(len(eods)): # for schleife erstellt liste mit den indices, an denen im eod array eine 0 steht
        if eods[i] == 0:
            index_is_zero.append(i)



    eods = eods.tolist() #befehl wandelt numpy array in eine liste um
    temperatures = temperatures.tolist()
    dates = dates.tolist()


    eods = [i for j, i in enumerate(eods) if j not in index_is_zero] #delets elements of the eod list which are zero
    temperatures = [i for j, i in enumerate(temperatures) if j not in index_is_zero] #delets the temperatures that belong to the zero eods
    dates = [i for j, i in enumerate(dates) if j not in index_is_zero] #delets the dates that belong to the zero eods


    trial_eods = OrderedDict()
    for d, t in zip(dates,eods):  #times and dates get ordered in a dictionary. therefore the date is the key on that the times of the date can be accesed
        if d not in trial_eods.keys():
            trial_eods[d] = []
        trial_eods[d].append(t)
    ticks = range(len(trial_eods.keys()))  #ticks are the scale of the x axis
    mean_eods = []
    std_eods = []
    for k in trial_eods.keys():
        mean_eods.append(np.mean(trial_eods[k]))
        std_eods.append(np.std(trial_eods[k]))  #np.std is the command for standard deviation

    fig = plt.figure()
    ax1 = fig.add_subplot(111)  # 111 means: one  x axis, one y axis, in one place. the tool subplot enables to put several diagramms in one plot
    ax1.scatter(ticks,mean_eods, color='g')  # ax.scatter is more or less the same like plt.scatter and creates a scatter plot (ax is the 'objektbezogene variante' while plt.scatter is originally taken from the matlab method)
    plt.errorbar(ticks, mean_eods,std_eods, color='g')  #plt.errorbar puts errorbars into the plot. plt.errorbar(x,y,standart deviation)
    ax1.set_xticks(ticks)
    xtickNames = ax1.set_xticklabels(trial_eods.keys())
    plt.setp(xtickNames, rotation=90, ha='center', fontsize=10) #rotates x labels for 45 degrees
    plt.subplots_adjust(left=0.1, bottom=0.15, right=0.9)
    ax1.set_ylabel('EOD [Hz]', color='g')
    ax1.tick_params(axis='y', which='major', direction ='inout', length=6, width=2)
    for tl in ax1.get_yticklabels():
        tl.set_color('g')

    trial_temperatures = OrderedDict()
    for d, t in zip(dates, temperatures):  #times and dates get ordered in a dictionary. therefore the date is the key on that the times of the date can be accesed
        if d not in trial_temperatures.keys():
            trial_temperatures[d] = []
        trial_temperatures[d].append(t)
    mean_temperatures = []
    std_temperatures = []
    for k in trial_temperatures.keys():
        mean_temperatures.append(np.mean(trial_temperatures[k]))
        std_temperatures.append(np.std(trial_temperatures[k]))  #np.std is the command for standard deviation
    new_style = {'grid': False}
    matplotlib.rc('axes', **new_style)
    for ticks in ax1.xaxis.get_ticklines() + ax1.yaxis.get_ticklines():
        ticks.set_color('green')

    ax2 = ax1.twinx() #creats second y-axis
    ticks2 = range(len(trial_temperatures.keys()))  #ticks are the scale of the x axis
    for ticks in ax2.xaxis.get_ticklines() + ax2.yaxis.get_ticklines():
        ticks.set_color('red')
    ax2.tick_params(axis='y', which='major', direction ='inout', length=6, width=2)
    ax2.scatter(ticks2,mean_temperatures, color='r')  # ax.scatter is more or less the same like plt.scatter and creates a scatter plot (ax is the 'objektbezogene variante' while plt.scatter is originally taken from the matlab method)
    plt.errorbar(ticks2, mean_temperatures,std_temperatures, color='r')  #plt.errorbar puts errorbars into the plot. plt.errorbar(x,y,standart deviation)
    ax2.set_ylabel(u'Temperatur [\xb0 C]', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')
    fig.set_size_inches(12, 6, forward=True)
    plt.savefig('date_temperature_eod_plot' + fish + '.pdf')
    plt.close()



def mean_eod_temperature_plot(eods, temperatures, fish):

    index_is_zero = []
    for i in np.arange(len(eods)): # for schleife erstellt liste mit den indices, an denen im eod array eine 0 steht
        if eods[i] == 0:
            index_is_zero.append(i)

    eods = eods.tolist() #befehl wandelt numpy array in eine liste um
    temperatures = temperatures.tolist()

    eods = [i for j, i in enumerate(eods) if j not in index_is_zero] #delets elements of the eod list which are zero
    temperatures = [i for j, i in enumerate(temperatures) if j not in index_is_zero] #delets the temperatures that belong to the zero eods

    trial_eods = dict()
    for d, t in zip(temperatures,eods):  #times and dates get ordered in a dictionary. therefore the date is the key on that the times of the date can be accesed
        if d not in trial_eods.keys():
            trial_eods[d] = []
        trial_eods[d].append(t)
    ticks = range(len(trial_eods.keys()))  #ticks are the scale of the x axis
    mean_eods = []
    std_eods = []
    for k in sorted(trial_eods.iterkeys()):
        mean_eods.append(np.mean(trial_eods[k]))
        std_eods.append(np.std(trial_eods[k]))  #np.std is the command for standard deviation

    fig = plt.figure()
    ax1 = fig.add_subplot(111)  # 111 means: one  x axis, one y axis, in one place. the tool subplot enables to put several diagramms in one plot
    ax1.scatter(ticks,mean_eods)  # ax.scatter is more or less the same like plt.scatter and creates a scatter plot (ax is the 'objektbezogene variante' while plt.scatter is originally taken from the matlab method)
    plt.errorbar(ticks, mean_eods,std_eods)  #plt.errorbar puts errorbars into the plot. plt.errorbar(x,y,standart deviation)
    ax1.set_xticks(ticks)
    xtickNames = ax1.set_xticklabels(sorted(trial_eods.iterkeys()))
    plt.setp(xtickNames, rotation=90, fontsize=10) #rotates x labels for 45 degrees
    plt.subplots_adjust(left=0.1, bottom=0.15, right=0.9)
    ax1.set_ylabel('EOD [Hz]')
    ax1.set_xlabel(u'Temperatur [\xb0 C]')
    #plt.title('Temperatur und EOD ' + fish)
    plt.savefig('temperature_mean_eod_plot' + fish + '.pdf')
    plt.close()

def eod_boxplot(eod1, eod2, eod3, eod4, eod5, eod6,fish1, fish2, fish3, fish4, fish5, fish6 ):

    fishnames = ['Fisch 1', 'Fisch 2', 'Fisch 3', 'Fisch 4', 'Fisch 5', 'Fisch 6']
    index_is_zero = []
    for i in np.arange(len(eod3)): # for schleife erstellt liste mit den indices, an denen im eod array eine 0 steht
        if eod3[i] == 0:
            index_is_zero.append(i)



    eod3 = eod3.tolist() #befehl wandelt numpy array in eine liste um
    eod3 = [i for j, i in enumerate(eod3) if j not in index_is_zero] #delets elements of the eod list which are zero



    eod = []
    eod.append(eod1)
    eod.append(eod2)
    eod.append(eod3)
    eod.append(eod4)
    eod.append(eod5)
    eod.append(eod6)




    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    boxplot_dict = ax.boxplot(eod)

    for b in boxplot_dict['boxes']:
        b.set_linewidth(1.5)
    for b in boxplot_dict['medians']:
        b.set_linewidth(1.5)

    ax.set_ylim((800, 1300))
    ax.set_xticklabels(fishnames, rotation=90)
    fig.subplots_adjust(bottom=0.3, wspace=0.3)
    ax.set_ylabel('EOD [Hz]')
    ax.yaxis.set_ticks(np.arange(800, 1350, 50))

    fig.savefig('eod_boxplot.pdf')
    plt.close()

    '''
    print np.median(eod1)
    print np.median(eod2)
    print np.median(eod3)
    print np.median(eod4)
    print np.median(eod5)
    print np.median(eod6)
    '''


def successrate_trials(successful, unsuccessful, fish, trial_number):

    g = float(len(successful)) # base value (Grundwert/100%) is the sum of all trials (successful + unsuccessful)

    percentage_successfull_trials = float(sum(successful))/g *100

    percentage_unsuccessfull_trials = float(sum(unsuccessful))/g *100

    return percentage_successfull_trials, percentage_unsuccessfull_trials, g



def succesrate_barplot(succesfull_percentage, unsuccesfull_percentage, trials_count, fish):

    print trials_count


    N = len(succesfull_percentage)
    ind = np.arange(N)  # the x locations for the groups
    width = 0.5  # the width of the bars
    fig, ax = plt.subplots()

    ## the bars
    first_bar = ax.bar(ind, succesfull_percentage, width, align='center')  # creates successful bar in green
    #second_bar = ax.bar(ind, unsuccesfull_percentage, width, color='red', align='center')

    labels = [u'Eingewöhnung', u'Konditioniereung auf S+', u'Gewöhnung an S-']

    def autolabel(bar, n):
        for i in np.arange(len(bar)):
            height = bar[i].get_height()
            ax.text(bar[i].get_x()+bar[i].get_width()/2., 1.01*height, '%d'%n[i], ha='center', va='bottom', fontsize=10)


    if len(trials_count) > 0:
        autolabel(first_bar, trials_count)

    ax.set_ylabel(u'erfolgreiche Versuchsdurchläufe [%]')
    ax.set_xticks(np.arange(len(succesfull_percentage)))
    ax.set_xticklabels(labels, ha='center')
    plt.ylim(0, 105)
        #p1 = "{:2.4}".format(str(p1))
        #p2 = "{:2.4}".format(str(p2))
        #ax.text(0.04, 50, p1 + ' %')
        #ax.text(0.14, 50, p2 + ' %')
    plt.savefig('Erfolgsquote' + fish + '.pdf')
    plt.show()


def date_mean_time_plot(dates, times, fish, trial_number):
    """
    function creates a plot per fish (4 fish) that shows the trial date on the x-axis and the mean time that the fish needed to complete on the y-axis
    :param dates: dates of successful trials
    :param times: times of successful trials, that means the time between opening the startbox and eating the larve
    :param fish: contains fish id
    :param trial_number: number of trial (versuch1-4)
    :return: plot that shows the mean times that the fish needed when the trial was successful for every trial date
    """

    if trial_number == 4 or trial_number == []:
        return
    else:
        trial_number = trial_number[0]
    trial_times = OrderedDict()
    for d, t in zip(dates,
                    times):  #times and dates get ordered in a dictionary. therefore the date is the key on that the times of the date can be accesed
        if d not in trial_times.keys():
            trial_times[d] = []
        trial_times[d].append(t)
    ticks = range(len(trial_times.keys()))  #ticks are the scale of the x axis
    mean_times = []
    std_times = []
    for k in trial_times.keys():
        mean_times.append(np.mean(trial_times[k]))
        std_times.append(np.std(trial_times[k]))  #np.std is the command for standard deviation

    fig = plt.figure()
    ax = fig.add_subplot(
        111)  # 111 means: one  x axis, one y axis, in one place. the tool subplot enables to put several diagramms in one plot
    ax.autoscale(enable=False, axis="x")
    ax.scatter(ticks,
               mean_times)  # ax.scatter is more or less the same like plt.scatter and creates a scatter plot (ax is the 'objektbezogene variante' while plt.scatter is originally taken from the matlab method)
    plt.errorbar(ticks, mean_times,
                 std_times)  #plt.errorbar puts errorbars into the plot. plt.errorbar(x,y,standart deviation)
    ax.set_xticks(ticks)
    xtickNames = ax.set_xticklabels(trial_times.keys())
    plt.setp(xtickNames, rotation=90, ha='center', fontsize=10) #rotates x labels for 45 degrees
    plt.subplots_adjust(left=0.1, bottom=0.15, right=0.95)
    plt.ylabel('Zeit [min]')
    plt.ylim(0, 25)
    plt.savefig('date_mean_time_plot_versuch' + trial_number + fish + '.pdf')
    plt.close()


def temperature_eod_plot(temperature, eod, fish):
    """
    function creates a scattter plot that shows the relation between eod and temperature. Furthermore a regression line is included
    :param temperature: array of temperatures of every trial
    :param eod: array of eods
    :param fish: fish id
    :return: plot with regression line
    """

    index_is_zero = []

    for i in np.arange(len(eod)): # for schleife erstellt liste mit den indices, an denen im eod array eine 0 steht (sinn der sache: falsche messdaten loswerden, wo eod null ist)
        if eod[i] == 0:
            index_is_zero.append(i)

    eod = eod.tolist() #befehl wandelt numpy array in eine liste um
    temperature = temperature.tolist()

    eod = [i for j, i in enumerate(eod) if j not in index_is_zero] #delets elements of the eod list which are zero
    temperature = [i for j, i in enumerate(temperature) if j not in index_is_zero] #delets the temperatures that belong to the zero eods

    fig = plt.figure()
    ax = fig.add_subplot(111)

    plt.xlabel(u'Temperatur [\xb0 C]')
    plt.ylabel('EOD [Hz]')
    m,y_achsenabschnitt, r_value, p_value, std_err = linregress(temperature,eod) #linregress gives steigung und y_achsenabschnitt der regressionsgeraden,pearsons r, p-wert and standartabweichung

    ax.plot(temperature,eod, 'yo', color = 'darkblue', label='Pearson\'r = %.3f, p = %.2e'%(r_value, p_value)) #creates a scatterplot with a label that includes pearsons r and p-value
    y =[] #the following lines create the regression line
    for t in temperature:
        f = m * t + y_achsenabschnitt
        y.append(f)
    plt.plot(temperature, y, color = 'red') #plots regression line
    plt.legend(loc=2, numpoints=1, markerscale=0., frameon=False)
    plt.savefig('eod_temperatur_plot ' + fish + '.pdf')
    plt.close()


def analyse_data(fish, file):
    """
    funktion gets the data that is already ordered after fish. Therefore the function gets four different datasets of four different fish
    this function should order the data into the different categories. therefore it goes through the data and creates lists

    :param fish: list that contains all data from one fish
    :return: orders the dates for every fish into different lists: time, eod, temperature, date, conductivity, successful, unsuccessful, time_successful, time_unsuccessful, date_successful, date_unsuccessful
    """

    time = []
    time_successful = []
    time_unsuccessful = []
    eod = []
    temperature = []
    date = []
    date_successful = []
    date_unsuccessful =[]
    conductivity = []
    successful = []
    unsuccessful = []
    trial_number = []


    for f in fish:
        parts = f.strip().split(',')  #f.split(',') divides the data that is currently organised in lines at the comma
        # f.strip takes away the /n on every end of a line

        if file == '/home/plappert/Documents/Versuchsprotokoll_Versuch4.csv':
            temperature. append(float(parts[6]))
            eod.append(float(parts[4]))
            date.append(parts[0])
            trial_number.append(parts[1])

        else:

            date.append(parts[0])  # command appends the first element of the line to the list date
            time.append(float(parts[3]))  #float converts strings to numbers
            eod.append(float(parts[4]))
            conductivity.append(float(parts[5]))
            temperature.append(float(parts[6]))
            trial_number.append(parts[1])
            if 'j' in parts[7]:
                successful.append(1)
                unsuccessful.append(0)
                time_successful.append(float(parts[3]))
                date_successful.append(parts[0])
            if 'n' in parts[7]:
                unsuccessful.append(1)
                successful.append(0)
                time_unsuccessful.append(float(parts[3]))
                date_unsuccessful.append(parts[0])

    return time, eod, temperature, date, conductivity, successful, unsuccessful, time_successful, time_unsuccessful, date_successful, date_unsuccessful, trial_number


def read_data(filename):
    """

    :param filename: contains the filename that should be read
    :return: returns one list per fish that contains all data of one fish
    """
    f = file(filename)
    lines = f.readlines()  # reads csv file and gives it back as lines, the columns of the file are divided through comma
    f.close()  # file hast to be closed after reading it
    chip = []  # empty list named after the first fish
    chap = []
    alfons = []
    trixi = []
    krummschwanz = []
    hermes = []

    for l in lines:
        if '2015albi02' in l:  # all lines that ingredient 2015albi02 get appended to the list chap
            chip.append(l)
        if '2015albi01' in l:
            chap.append(l)
        if '2014albi08' in l:
            alfons.append(l)
        if '2013albi14' in l:
            trixi.append(l)
        if '2013albi09' in l:
            krummschwanz.append(l)
        if '2012albi01' in l:
            hermes.append(l)

    return chip, chap, alfons, trixi, krummschwanz, hermes


if __name__ == '__main__':  # the code doesn't run if someone is adding it to his or her code it only runs if it is the main code
    files=glob.glob('/home/plappert/Documents/Versuchsprotokoll*.csv')
    print files


    fish1 = '2015albi02'
    fish2 = '2015albi01'
    fish3 = '2014albi08'
    fish4 = '2013albi14'
    fish5 = '2013albi09'
    fish6 = '2012albi01'


    temperature_chip = [np.array([]) for e in np.arange(len(files))]
    eod_chip = [np.array([]) for e in np.arange(len(files))]
    temperature_chap = [np.array([]) for e in np.arange(len(files))]
    eod_chap = [np.array([]) for e in np.arange(len(files))]
    temperature_alfons= [np.array([]) for e in np.arange(len(files))]
    eod_alfons = [np.array([]) for e in np.arange(len(files))]
    temperature_trixi = [np.array([]) for e in np.arange(len(files))]
    eod_trixi = [np.array([]) for e in np.arange(len(files))]
    temperature_krummschwanz = [np.array([]) for e in np.arange(len(files))]
    eod_krummschwanz = [np.array([]) for e in np.arange(len(files))]
    temperature_hermes = [np.array([]) for e in np.arange(len(files))]
    eod_hermes = [np.array([]) for e in np.arange(len(files))]


    date_chip = [np.array([]) for e in np.arange(len(files))]
    date_chap = [np.array([]) for e in np.arange(len(files))]
    date_alfons = [np.array([]) for e in np.arange(len(files))]
    date_trixi = [np.array([]) for e in np.arange(len(files))]
    date_krummschwanz = [np.array([]) for e in np.arange(len(files))]
    date_hermes = [np.array([]) for e in np.arange(len(files))]

    all_trials_succesfull_percentage1 = []
    all_trials_unsuccesfull_percentage1 = []
    all_trials_succesfull_percentage2 = []
    all_trials_unsuccesfull_percentage2 = []
    all_trials_succesfull_percentage3 = []
    all_trials_unsuccesfull_percentage3 = []
    all_trials_succesfull_percentage4 = []
    all_trials_unsuccesfull_percentage4 = []

    all_trials_count1 = []
    all_trials_count2 = []
    all_trials_count3 = []
    all_trials_count4 = []

    for enu, f in enumerate(files):
        print f
        chip, chap, alfons, trixi, krummschwanz, hermes = read_data(f)  # Funktion, die die Daten fuer jeden der vier Fische auslesen soll

        sns.despine()

        time1, EOD1, temperature1, date1, conductiyity1, successful1, unsuccessful1, time_successful1, time_unsuccesful1, date_successful1, date_unsuccessful1, trial_number1 = analyse_data(chip,f)  #the data of the function read_data, that already ordered the data after fish, get devided into different lists like date, time etc
        time2, EOD2, temperature2, date2, conductivity2, successful2, unsuccessful2, time_successful2, time_unsuccesful2, date_successful2, date_unsuccessful2, trial_number2 = analyse_data(chap,f)
        time3, EOD3, temperature3, date3, conductivity3, successful3, unsuccessful3, time_successful3, time_unsuccesful3, date_successful3, date_unsuccessful3, trial_number3 = analyse_data(alfons,f)
        time4, EOD4, temperature4, date4, conductivity4, successful4, unsuccessful4, time_successful4, time_unsuccesful4, date_successful4, date_unsuccessful4, trial_number4 = analyse_data(trixi,f)
        time5, EOD5, temperature5, date5, conductivity5, successful5, unsuccessful5, time_successful5, time_unsuccesful5, date_successful5, date_unsuccessful5, trial_number5 = analyse_data(krummschwanz,f)
        time6, EOD6, temperature6, date6, conductivity6, successful6, unsuccessful6, time_successful6, time_unsuccesful6, date_successful6, date_unsuccessful6, trial_number6 = analyse_data(hermes,f)


        if trial_number1 == []:
            number = float(trial_number6[0]) - 1
            number = int(number)
        else:
            number = float(trial_number1[0]) - 1
            number = int(number)

        temperature_chip[number] = temperature1
        eod_chip[number] = EOD1
        date_chip[number] = date1

        temperature_chap[number] = temperature2
        eod_chap[number] = EOD2
        date_chap[number] = date2


        temperature_alfons[number] = temperature3
        eod_alfons[number] = EOD3
        date_alfons[number] = date3


        temperature_trixi[number] = temperature4
        eod_trixi[number] = EOD4
        date_trixi[number] = date4

        temperature_krummschwanz[number] = temperature5
        eod_krummschwanz[number] = EOD5
        date_krummschwanz[number] = date5

        temperature_hermes[number] = temperature6
        eod_hermes[number] = EOD6
        date_hermes[number] = date6

        if number < 3:
            percentage_succesful1, percentage_unsuccesful1, trials_count1 = successrate_trials(successful1, unsuccessful1, fish1, trial_number1)
            percentage_succesful2, percentage_unsuccesful2, trials_count2 = successrate_trials(successful2, unsuccessful2, fish2, trial_number2)
            percentage_succesful3, percentage_unsuccesful3, trials_count3 = successrate_trials(successful3, unsuccessful3, fish3, trial_number3)
            percentage_succesful4, percentage_unsuccesful4, trials_count4 = successrate_trials(successful4, unsuccessful4, fish4, trial_number4)

            all_trials_succesfull_percentage1.append(percentage_succesful1)
            all_trials_unsuccesfull_percentage1.append(percentage_unsuccesful1)
            all_trials_succesfull_percentage2.append(percentage_succesful2)
            all_trials_unsuccesfull_percentage2.append(percentage_unsuccesful2)
            all_trials_succesfull_percentage3.append(percentage_succesful3)
            all_trials_unsuccesfull_percentage3.append(percentage_unsuccesful3)
            all_trials_succesfull_percentage4.append(percentage_succesful4)
            all_trials_unsuccesfull_percentage4.append(percentage_unsuccesful4)

            all_trials_count1.append(trials_count1)
            all_trials_count2.append(trials_count2)
            all_trials_count3.append(trials_count3)
            all_trials_count4.append(trials_count4)

        date_mean_time_plot(date_successful1, time_successful1, fish1, trial_number1)  #function plots the mean times the fish needed for the trials per day including errorbars, but only for successful trials
        date_mean_time_plot(date_successful2, time_successful2, fish2, trial_number2)
        date_mean_time_plot(date_successful3, time_successful3, fish3, trial_number3)
        date_mean_time_plot(date_successful4, time_successful4, fish4, trial_number4)
        date_mean_time_plot(date_successful6, time_successful6, fish6, trial_number6)



    succesrate_barplot(all_trials_succesfull_percentage1, all_trials_unsuccesfull_percentage1, all_trials_count1, fish1)
    succesrate_barplot(all_trials_succesfull_percentage2, all_trials_unsuccesfull_percentage2, all_trials_count2, fish2)
    succesrate_barplot(all_trials_succesfull_percentage3, all_trials_unsuccesfull_percentage3, all_trials_count3, fish3)
    succesrate_barplot(all_trials_succesfull_percentage4, all_trials_unsuccesfull_percentage4, all_trials_count4, fish4)

    temperature_chip = np.hstack(temperature_chip)
    temperature_chap = np.hstack(temperature_chap)
    temperature_alfons = np.hstack(temperature_alfons)
    temperature_trixi = np.hstack(temperature_trixi)
    temperature_krummschwanz = np.hstack(temperature_krummschwanz)
    temperature_hermes = np.hstack(temperature_hermes)

    eod_chip = np.hstack(eod_chip)
    eod_chap = np.hstack(eod_chap)
    eod_alfons = np.hstack(eod_alfons)
    eod_trixi = np.hstack(eod_trixi)
    eod_krummschwanz = np.hstack(eod_krummschwanz)
    eod_hermes = np.hstack(eod_hermes)

    date_chip = np.hstack(date_chip)
    date_chap = np.hstack(date_chap)
    date_alfons = np.hstack(date_alfons)
    date_trixi = np.hstack(date_trixi)
    date_krummschwanz = np.hstack(date_krummschwanz)
    date_hermes = np.hstack(date_hermes)


    temperature_eod_plot(temperature_chip, eod_chip, fish1)  #function plots eod(y-axis) in relation to temperature(x-axis)
    temperature_eod_plot(temperature_chap, eod_chap, fish2)
    temperature_eod_plot(temperature_alfons, eod_alfons, fish3)
    temperature_eod_plot(temperature_trixi, eod_trixi, fish4)
    temperature_eod_plot(temperature_krummschwanz, eod_krummschwanz, fish5)
    temperature_eod_plot(temperature_hermes, eod_hermes, fish6)


    mean_eod_temperature_plot(eod_chip, temperature_chip, fish1)
    mean_eod_temperature_plot(eod_chap, temperature_chap, fish2)
    mean_eod_temperature_plot(eod_alfons, temperature_alfons, fish3)
    mean_eod_temperature_plot(eod_trixi, temperature_trixi, fish4)
    mean_eod_temperature_plot(eod_krummschwanz, temperature_krummschwanz, fish5)
    mean_eod_temperature_plot(eod_hermes, temperature_hermes, fish6)

    eod_boxplot(eod_chip, eod_chap, eod_alfons, eod_trixi, eod_krummschwanz, eod_hermes, fish1, fish2, fish3, fish4, fish5, fish6)

    date_eod_temperature_plot(date_chip, eod_chip, temperature_chip, fish1)
    date_eod_temperature_plot(date_chap, eod_chap, temperature_chap, fish2)
    date_eod_temperature_plot(date_alfons, eod_alfons, temperature_alfons, fish3)
    date_eod_temperature_plot(date_trixi, eod_trixi, temperature_trixi, fish4)
    date_eod_temperature_plot(date_krummschwanz, eod_krummschwanz, temperature_krummschwanz, fish5)
    date_eod_temperature_plot(date_hermes, eod_hermes, temperature_hermes, fish6)






