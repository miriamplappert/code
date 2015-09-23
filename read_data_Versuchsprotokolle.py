#! /usr/bin/env python
__author__ = 'plappert'
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from IPython import embed
import glob
from scipy.stats import linregress


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
    ax1.scatter(ticks,mean_eods, color='b')  # ax.scatter is more or less the same like plt.scatter and creates a scatter plot (ax is the 'objektbezogene variante' while plt.scatter is originally taken from the matlab method)
    plt.errorbar(ticks, mean_eods,std_eods, color='b')  #plt.errorbar puts errorbars into the plot. plt.errorbar(x,y,standart deviation)
    ax1.set_xticks(ticks)
    xtickNames = ax1.set_xticklabels(trial_eods.keys())
    plt.setp(xtickNames, rotation=45, fontsize=10) #rotates x labels for 45 degrees
    plt.subplots_adjust(left=0.1, bottom=0.15, right=0.9)
    ax1.set_ylabel('EOD [Millivolt]', color='b')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')


    trial_temperatures = OrderedDict()
    for d, t in zip(dates, temperatures):  #times and dates get ordered in a dictionary. therefore the date is the key on that the times of the date can be accesed
        if d not in trial_temperatures.keys():
            trial_temperatures[d] = []
        trial_temperatures[d].append(t)
    ticks = range(len(trial_temperatures.keys()))  #ticks are the scale of the x axis
    mean_temperatures = []
    std_temperatures = []
    for k in trial_temperatures.keys():
        mean_temperatures.append(np.mean(trial_temperatures[k]))
        std_temperatures.append(np.std(trial_temperatures[k]))  #np.std is the command for standard deviation

    ax2 = ax1.twinx() #creats second y-axis
    ax2.scatter(ticks,mean_temperatures, color='g')  # ax.scatter is more or less the same like plt.scatter and creates a scatter plot (ax is the 'objektbezogene variante' while plt.scatter is originally taken from the matlab method)
    plt.errorbar(ticks, mean_temperatures,std_temperatures, color='g')  #plt.errorbar puts errorbars into the plot. plt.errorbar(x,y,standart deviation)
    ax2.set_ylabel('Temperatur [Grad Celsius]', color='g')
    for tl in ax2.get_yticklabels():
        tl.set_color('g')
    ax1.set_axis_bgcolor('powderblue')
    plt.grid(color='white', linestyle='-')
    ax1.set_axisbelow(True)
    plt.xlabel('Datum')
    plt.title('Korrelation von Temperatur und EOD ' + fish)
    plt.savefig('date_temperature_eod_plot' + fish + '.pdf')
    plt.show()


def successrate_bar_plot(successful, unsuccessful, fish, trial_number, dates):
    """
    function shows how many trials in every step (versuch1-4) were successful and how many were unsuccessful.
    therefor the function creates for every of the four fish a bar plot
    :param fish: contains the fish-id (eg.: 2015albi01)
    """
    print trial_number
    trial_number = trial_number[0]
    g = len(successful) # base value (Grundwert/100%) is the sum of all trials (successful + unsuccessful)
    if g == 0:
        return
    else:

        successful_trials_dates = OrderedDict() #ordered dicitonary that contains for every date the amount of zeros and ones of the right choice list
        for d, t in zip(dates,successful):
            if d not in successful_trials_dates.keys():
                successful_trials_dates[d] = []
            successful_trials_dates[d].append(t)


        unsuccessful_trials_dates = OrderedDict()
        for e, u in zip(dates,unsuccessful):  #times and dates get ordered in a dictionary. therefore the date is the key on that the times of the date can be accesed
            if e not in unsuccessful_trials_dates.keys():
                unsuccessful_trials_dates[e] = []
            unsuccessful_trials_dates[e].append(u)

        print successful_trials_dates
        print unsuccessful_trials_dates

        percentage_successfull_trials = []
        percentage_unsuccessfull_trials = []

        for z in successful_trials_dates.keys(): #for loop goes through the dicitonary keys
            x = successful_trials_dates[z]
            y = float(sum(x)) / float(len(x)) *100
            percentage_successfull_trials.append(y)

        for a in unsuccessful_trials_dates.keys():
            b = unsuccessful_trials_dates[a]
            c = float(sum(b)) / float(len(b)) *100
            percentage_unsuccessfull_trials.append(c)

        N = len(unsuccessful_trials_dates.keys())
        ind = np.arange(N)  # the x locations for the groups
        width = 0.5  # the width of the bars
        fig, ax = plt.subplots()

        ## the bars
        first_bar = ax.bar(ind, percentage_successfull_trials, width, color='blue')  # creates successful bar in green

        ax.set_ylabel('Erfolgreiche Trials [%]')
        ax.set_title('Erfolgsquote Versuch ' + trial_number + ' ' + fish)
        ax.set_xticks(ind + width)
        ax.set_xticklabels('Versuch ' + trial_number)
        plt.ylim(0, 100)
        #p1 = "{:2.4}".format(str(p1))
        #p2 = "{:2.4}".format(str(p2))
        #ax.text(0.04, 50, p1 + ' %')
        #ax.text(0.14, 50, p2 + ' %')
        plt.setp(ax.get_xticklabels(), visible=False) # lets x tick labels disappear
        ax.set_axis_bgcolor('powderblue')
        plt.grid(color='white', linestyle='-')
        ax.set_axisbelow(True)

        plt.savefig('Erfolgsquote_Versuch ' + trial_number + fish + '.pdf')
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
    ax.scatter(ticks,
               mean_times)  # ax.scatter is more or less the same like plt.scatter and creates a scatter plot (ax is the 'objektbezogene variante' while plt.scatter is originally taken from the matlab method)
    plt.errorbar(ticks, mean_times,
                 std_times)  #plt.errorbar puts errorbars into the plot. plt.errorbar(x,y,standart deviation)
    ax.set_xticks(ticks)
    xtickNames = ax.set_xticklabels(trial_times.keys())
    plt.setp(xtickNames, rotation=45, fontsize=10) #rotates x labels for 45 degrees
    plt.subplots_adjust(left=0.1, bottom=0.15, right=0.95)
    plt.xlabel('Datum')
    plt.ylabel('Zeit [Minuten]')
    plt.title('Mittlere Zeit Versuch ' +trial_number + ' ' + fish)
    plt.ylim(0, 25)
    ax.set_axis_bgcolor('powderblue')
    plt.grid(color='white', linestyle='-')
    ax.set_axisbelow(True)

    plt.savefig('date_mean_time_plot_versuch' + trial_number + fish + '.pdf')
    plt.show()


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

    plt.xlabel('Temperatur [Grad Celsius]')
    plt.ylabel('EOD [Millivolt]')
    plt.title('Relation von Temperatur und EOD ' + fish)
    m,y_achsenabschnitt, r_value, p_value, std_err = linregress(temperature,eod,) #linregress gives steigung und y_achsenabschnitt der regressionsgeraden,pearsons r, p-wert and standartabweichung

    ax.plot(temperature,eod, 'yo', color='0.75', label='Pearson\'r = %.3f, p = %.2e'%(r_value, p_value)) #creates a scatterplot with a label that includes pearsons r and p-value
    y =[] #the following lines create the regression line
    for t in temperature:
        f = m * t + y_achsenabschnitt
        y.append(f)
    plt.plot(temperature, y) #plots regression line
    ax.set_axis_bgcolor('powderblue')
    plt.grid(color='white', linestyle='-')
    ax.set_axisbelow(True)
    plt.legend(loc=2, numpoints=1, markerscale=0., frameon=False)
    plt.savefig('eod_temperatur_plot ' + fish + '.pdf')
    plt.show()


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
    for l in lines:
        if '2015albi02' in l:  # all lines that ingredient 2015albi02 get appended to the list chap
            chip.append(l)
        if '2015albi01' in l:
            chap.append(l)
        if '2014albi08' in l:
            alfons.append(l)
        if '2013albi14' in l:
            trixi.append(l)
    return chip, chap, alfons, trixi


if __name__ == '__main__':  # the code doesn't run if someone is adding it to his or her code it only runs if it is the main code
    files=glob.glob('/home/plappert/Documents/Versuchsprotokoll*.csv')
    print files

    temperature_chip = [np.array([]) for e in np.arange(len(files))]
    eod_chip = [np.array([]) for e in np.arange(len(files))]
    temperature_chap = [np.array([]) for e in np.arange(len(files))]
    eod_chap = [np.array([]) for e in np.arange(len(files))]
    temperature_alfons= [np.array([]) for e in np.arange(len(files))]
    eod_alfons = [np.array([]) for e in np.arange(len(files))]
    temperature_trixi = [np.array([]) for e in np.arange(len(files))]
    eod_trixi = [np.array([]) for e in np.arange(len(files))]

    date_chip = [np.array([]) for e in np.arange(len(files))]
    date_chap = [np.array([]) for e in np.arange(len(files))]
    date_alfons = [np.array([]) for e in np.arange(len(files))]
    date_trixi = [np.array([]) for e in np.arange(len(files))]

    for enu, f in enumerate(files):
        chip, chap, alfons, trixi = read_data(f)  # Funktion, die die Daten fuer jeden der vier Fische auslesen soll
        time1, EOD1, temperature1, date1, conductiyity1, successful1, unsuccessful1, time_successful1, time_unsuccesful1, date_successful1, date_unsuccessful1, trial_number1 = analyse_data(chip,f)  #the data of the function read_data, that already ordered the data after fish, get devided into different lists like date, time etc
        time2, EOD2, temperature2, date2, conductivity2, successful2, unsuccessful2, time_successful2, time_unsuccesful2, date_successful2, date_unsuccessful2, trial_number2 = analyse_data(chap,f)
        time3, EOD3, temperature3, date3, conductivity3, successful3, unsuccessful3, time_successful3, time_unsuccesful3, date_successful3, date_unsuccessful3, trial_number3 = analyse_data(alfons,f)
        time4, EOD4, temperature4, date4, conductivity4, successful4, unsuccessful4, time_successful4, time_unsuccesful4, date_successful4, date_unsuccessful4, trial_number4 = analyse_data(trixi,f)



        temperature_chip[enu] = temperature1
        temperature_chap[enu] = temperature2
        temperature_alfons[enu] = temperature3
        temperature_trixi[enu] = temperature4

        eod_chip[enu] = EOD1
        eod_chap[enu] = EOD2
        eod_alfons[enu] = EOD3
        eod_trixi[enu] = EOD4

        date_chip[enu] = date1
        date_chap[enu] = date2
        date_alfons[enu] = date3
        date_trixi[enu] = date4

        fish1 = '2015albi02'
        fish2 = '2015albi01'
        fish3 = '2014albi08'
        fish4 = '2013albi14'

        print trial_number3

        '''
        successrate_bar_plot(successful1, unsuccessful1, fish1, trial_number1, date1)  #function creates bar plots that show how many of the trials where successful
        successrate_bar_plot(successful2, unsuccessful2, fish2, trial_number2, date2)
        #successrate_bar_plot(successful3, unsuccessful3, fish3, trial_number3, date3)
        successrate_bar_plot(successful4, unsuccessful4, fish4, trial_number4, date4)

        date_mean_time_plot(date_successful1, time_successful1, fish1, trial_number1)  #function plots the mean times the fish needed for the trials per day including errorbars, but only for successful trials
        date_mean_time_plot(date_successful2, time_successful2, fish2, trial_number2)
        #date_mean_time_plot(date_successful3, time_successful3, fish3, trial_number3)
        date_mean_time_plot(date_successful4, time_successful4, fish4, trial_number4)
        '''

    temperature_chip = np.hstack(temperature_chip)
    temperature_chap = np.hstack(temperature_chap)
    temperature_alfons = np.hstack(temperature_alfons)
    temperature_trixi = np.hstack(temperature_trixi)

    eod_chip = np.hstack(eod_chip)
    eod_chap = np.hstack(eod_chap)
    eod_alfons = np.hstack(eod_alfons)
    eod_trixi = np.hstack(eod_trixi)

    date_chip = np.hstack(date_chip)
    date_chap = np.hstack(date_chap)
    date_alfons = np.hstack(date_alfons)
    date_trixi = np.hstack(date_trixi)

    temperature_eod_plot(temperature_chip, eod_chip, fish1)  #function plots eod(y-axis) in relation to temperature(x-axis)
    temperature_eod_plot(temperature_chap, eod_chap, fish2)
    temperature_eod_plot(temperature_alfons, eod_alfons, fish3)
    temperature_eod_plot(temperature_trixi, eod_trixi, fish4)

    date_eod_temperature_plot(date_chip, eod_chip, temperature_chip, fish1)
    date_eod_temperature_plot(date_chap, eod_chap, temperature_chap, fish2)
    date_eod_temperature_plot(date_alfons, eod_alfons, temperature_alfons, fish3)
    date_eod_temperature_plot(date_trixi, eod_trixi, temperature_trixi, fish4)












