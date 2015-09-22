__author__ = 'plappert'
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from IPython import embed
import scipy.stats




def binomialtest(right_choice, wrong_choice, no_choice):
    # Ha = der Fisch ebntscheidet sich oefter fuer die richtige elektrode (k1) als fuer die falsche elektrode (k2) --> Ha: k1 > k2
    # H0: k1 <= k2

    index_is_one = []

    for i in np.arange(len(no_choice)):
        if no_choice[i] == 1:
            index_is_one.append(i)

    right_choice = [i for j, i in enumerate(right_choice) if j not in index_is_one]
    wrong_choice = [i for j, i in enumerate(wrong_choice) if j not in index_is_one]

    k1 = sum(right_choice)
    n = len(wrong_choice) # n = number of trials
    p = 0.5 # p = The hypothesized probability of success. 0 <= p <= 1. The default value is p = 0.5



    p_value = scipy.stats.binom_test(k1,n,p) # k1 = number of right choices, n = number of trials, p = probability of right choices
    print p_value

    return p_value


def load(filename):
    """
    FUNKTION VON EILEEN

    Loads a data file saved by relacs. Returns a tuple of dictionaries
    containing the data and the header information

    :param filename: Filename of the data file.
    :type filename: string
    :returns:  a tuple of dictionaries containing the head information and the data.
    :rtype: tuple

    """
    with open(filename, 'r') as fid:
        L = [l.lstrip().rstrip() for l in fid.readlines()]

    ret = []
    dat = {}
    X = []
    keyon = False
    currkey = None
    for l in L:
        # if empty line and we have data recorded
        if (not l or l.startswith('#')) and len(X) > 0:
            keyon = False
            currkey = None
            dat['data'] = np.array(X)
            ret.append(dat)
            X = []
            dat = {}

        if '---' in l:
            continue
        if l.startswith('#'):
            if ":" in l:
                tmp = [e.rstrip().lstrip() for e in l[1:].split(':')]
                if currkey is None:
                    dat[tmp[0]] = tmp[1]
                else:
                    dat[currkey][tmp[0]] = tmp[1]
            elif "=" in l:
                tmp = [e.rstrip().lstrip() for e in l[1:].split('=')]
                if currkey is None:
                    dat[tmp[0]] = tmp[1]
                else:
                    dat[currkey][tmp[0]] = tmp[1]
            elif l[1:].lower().startswith('key'):
                dat['key'] = []

                keyon = True
            elif keyon:

                dat['key'].append(tuple([e.lstrip().rstrip() for e in l[1:].split()]))
            else:
                currkey = l[1:].rstrip().lstrip()
                dat[currkey] = {}

        elif l: # if l != ''
            keyon = False
            currkey = None
            X.append( [float(e) for e in l.split()])

    if len(X) > 0:
            dat['data'] = np.array(X)
    else:
            dat['data'] = []
    ret.append(dat)

    return tuple(ret)

def performance_without_no_choice(right_choice, wrong_choice, no_choice, dates, fish):

    index_is_one = []

    for i in np.arange(len(no_choice)): # for schleife erstellt liste mit den indices, an denen im eod array eine 0 steht (sinn der sache: falsche messdaten loswerden, wo eod null ist)
        if no_choice[i] == 1:
            index_is_one.append(i)

    right_choice = [i for j, i in enumerate(right_choice) if j not in index_is_one]
    wrong_choice = [i for j, i in enumerate(wrong_choice) if j not in index_is_one] #delets the temperatures that belong to the zero eods
    dates = [i for j, i in enumerate(dates) if j not in index_is_one]

    right_choices_dates = OrderedDict() #ordered dicitonary that contains for every date the amount of zeros and ones of the right choice list
    for d, t in zip(dates,right_choice):
        if d not in right_choices_dates.keys():
            right_choices_dates[d] = []
        right_choices_dates[d].append(t)


    wrong_choices_dates = OrderedDict()
    for e, u in zip(dates,wrong_choice):  #times and dates get ordered in a dictionary. therefore the date is the key on that the times of the date can be accesed
        if e not in wrong_choices_dates.keys():
            wrong_choices_dates[e] = []
        wrong_choices_dates[e].append(u)

    percentage_right_choices = []
    percentage_wrong_choices = []

    for z in right_choices_dates.keys(): #for loop goes through the dicitonary keys
        x = right_choices_dates[z]
        y = float(sum(x)) / float(len(x)) *100
        percentage_right_choices.append(y)

    for a in wrong_choices_dates.keys():
        b = wrong_choices_dates[a]
        c = float(sum(b)) / float(len(b)) *100
        percentage_wrong_choices.append(c)

    N = len(right_choices_dates.keys())
    ind = np.arange(N)  # the x locations for the groups
    width = 0.7  # the width of the bars
    fig, ax = plt.subplots()

    ## the bars
    first_bar = ax.bar(ind, percentage_right_choices, width, color='blue')

    ax.set_ylabel('Richtige Entscheidungen [%]')
    ax.set_title('Performance ' + fish)
    ax.set_xticks(ind)
    ax.set_xticklabels(right_choices_dates.keys(), rotation=45, fontsize=10)
    plt.axhline(y=50, xmin=0, xmax=1, hold=None, color='white', linewidth=2, linestyle='dashed')
    plt.ylim(0, 100)
    ax.set_axis_bgcolor('lightblue')
    plt.grid(color='white', linestyle='-')
    ax.set_axisbelow(True)
    plt.savefig('Performance_without_no_choice' + fish + '.pdf')
    plt.show()

    mean_performance = np.mean(percentage_right_choices)
    print mean_performance

    return mean_performance



def performance_line_plot(percentage_right_choices, dates, fish):
    '''

    :param percentage_right_choices: list of right choice percentages for every date
    :param dates: dates of trials
    :param fish: name of the fish
    :return: creates a line plot that shows what percentage of choices was right for every trial date
    '''
    ticks = range(len(dates))
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(ticks,percentage_right_choices)
    plt.plot(ticks, percentage_right_choices)
    plt.axhline(y=50, xmin=0, xmax=1, hold=None, color='black', linestyle='dashed')
    ax1.set_xticks(ticks)
    xtickNames = ax1.set_xticklabels(dates)
    plt.setp(xtickNames, rotation=45, fontsize=10) #rotates x labels for 45 degrees
    plt.subplots_adjust(left=0.1, bottom=0.15, right=0.9)
    ax1.set_ylabel('Performance [%]')
    plt.xlabel('Datum')
    plt.ylim(0,120)
    plt.title('Performance' + fish)
    plt.savefig('performance_line_plot' + fish + '.pdf')
    plt.show()


def performance_bar_plot(rewarded_electrodes, chosen_elektrodes, dates, fish):
    '''

    :param rewarded_electrodes: list of rewarded eectrodes of every fish for all trials
    :param chosen_elektrodes: list of the electrodes the fish chose of all trials
    :param dates: list of dates that fits to the rewarded and chosen electrodes list
    :param fish: name of the fish
    :return: creates a bar plot that shows in what percentage the fish chose the right electrode for every date, returns percentages of right, wrong and no choices and a list of the fitting dates
    '''

    right_choice = [] # empty lists that will be filled up with ones and zeros (1 means true, 0 means false)
    wrong_choice = []
    no_choice = []

    for i in np.arange(len(rewarded_electrodes)): # for loop goes through the indices of the rewarded_electrodes list
        if rewarded_electrodes[i] == chosen_elektrodes[i]: # right choice rewarded electrode = chosen electrode
            right_choice.append(1) # number 1 stands for true --> in this case the decision was right so I append a one to my right choice list
            wrong_choice.append(0) # zero stands for false
            no_choice.append(0)
        elif chosen_elektrodes[i] == 0: # in this case the fish didn't make a choice at all
            right_choice.append(0)
            wrong_choice.append(0)
            no_choice.append(1)
        else: # last case should be that the fish chose the wrong electrode (eg: rewarded electrode = 1, chosen electrode = 2)
            right_choice.append(0)
            wrong_choice.append(1)
            no_choice.append(0)

    right_choices_dates = OrderedDict() #ordered dicitonary that contains for every date the amount of zeros and ones of the right choice list
    for d, t in zip(dates,right_choice):
        if d not in right_choices_dates.keys():
            right_choices_dates[d] = []
        right_choices_dates[d].append(t)


    wrong_choices_dates = OrderedDict()
    for e, u in zip(dates,wrong_choice):  #times and dates get ordered in a dictionary. therefore the date is the key on that the times of the date can be accesed
        if e not in wrong_choices_dates.keys():
            wrong_choices_dates[e] = []
        wrong_choices_dates[e].append(u)


    no_choices_dates = OrderedDict()
    for g, v in zip(dates,no_choice):  #times and dates get ordered in a dictionary. therefore the date is the key on that the times of the date can be accesed
        if g not in no_choices_dates.keys():
            no_choices_dates[g] = []
        no_choices_dates[g].append(v)


    percentage_right_choices = []
    percentage_wrong_choices = []
    percentage_no_choices = []

    for z in right_choices_dates.keys(): #for loop goes through the dicitonary keys
        x = right_choices_dates[z]
        y = float(sum(x)) / float(len(x)) *100
        percentage_right_choices.append(y)

    for a in wrong_choices_dates.keys():
        b = wrong_choices_dates[a]
        c = float(sum(b)) / float(len(b)) *100
        percentage_wrong_choices.append(c)


    for j in no_choices_dates.keys():
        k = no_choices_dates[j]
        l = float(sum(k)) / float(len(k)) *100
        percentage_no_choices.append(l)



    N = len(right_choices_dates.keys())
    ind = np.arange(N)  # the x locations for the groups
    width = 0.1  # the width of the bars
    fig, ax = plt.subplots()

    ## the bars
    first_bar = ax.bar(ind, percentage_right_choices, width, color='green')  # creates successful bar in green
    second_bar = ax.bar(ind + width, percentage_wrong_choices, width, color='red')  # creates unsuccessful bar in red
    third_bar = ax.bar(ind + width, percentage_no_choices, width, color='grey')  # creates unsuccessful bar in red

    ax.set_ylabel('Prozensatz[%]')
    ax.set_title('Performance ' + fish)
    ax.set_xticks(ind+width)
    ax.set_xticklabels(right_choices_dates.keys(), rotation=45, fontsize=10)
    plt.axhline(y=50, xmin=0, xmax=1, hold=None, color='black', linestyle='dashed')
    plt.ylim(0, 120)
    ax.legend((first_bar[0], second_bar[0], third_bar[0]), ('richtige Elektrode', 'falsche Elektrode', 'keine Entscheidung'))
    plt.savefig('Performance' + fish + '.pdf')
    #plt.show()
    plt.close()

    return percentage_right_choices, percentage_wrong_choices, percentage_no_choices, right_choices_dates.keys(), right_choice, wrong_choice, no_choice


def read_relacs_files(relacs_files):
    '''
    this function reads the relacs file 'choice dat' and orders the data given in the file to lists.
    :param relacs_files: relacs files is a list of all relacs file names of one fish. some file names have additional this numbers like this _1, this shows in what line of the choices.dat file the data is
    :return: returns the following list: dicitonaries, keys, datas, rewarded_electrode, chosen_electrode, amplitude_1, amplitude_2, eod, electrode_frequency, start_time, choice_time
    '''

    trial_numbers = []
    dicitonaries = []
    keys = []
    datas = []

    signal_duration = []
    beat_frequency = []
    rewarded_electrode = []
    chosen_electrode = []
    amplitude_1 = []
    amplitude_2 = []
    eod = []
    electrode_frequency = []
    start_time = []
    choice_time = []

    for f in relacs_files:
        if len(f) > 13: #for filenames that include an addition like _1, _2 and so on
            trial_number = float(f[14:])
            trial_number += -1
            #trial number gibt an, der wievielte durchlauf es innerhalb der datei war, -1 damit indices stimmen
            f = f[0:13] # nur die ersten 13 zeichen sind der filename, _1, _2 usw. wird abgeschnitten


        else:
            trial_number = 0 #bei files die nur einen durchlauf hatten und somit keinen zusatz wie _1

        x = load('../../Data/relacs_files/' + f + '/choices.dat')  # x is a tuple; x[0] is the dictionary that I want.

        dictionary = x[0]
        keys.append(dictionary.keys())
        data = dictionary['data']
        if trial_number > 9:
            trial_number += -10
            dicitonary2 = x[1]
            data2 = dicitonary2['data']
            rewarded_electrode.append(data2[trial_number][0])
            chosen_electrode.append(data2[trial_number][1])
            amplitude_1.append(round(data2[trial_number][3],2))
            amplitude_2.append(round(data2[trial_number][4],2))
            eod.append(round(data2[trial_number][5]))
            electrode_frequency.append(round(data2[trial_number][6]))
            start_time.append(round(data2[trial_number][8]))
            choice_time.append(round(data2[trial_number][9]))
        else:
            rewarded_electrode.append(data[trial_number][0])
            chosen_electrode.append(data[trial_number][1])
            amplitude_1.append(round(data[trial_number][3],2))
            amplitude_2.append(round(data[trial_number][4],2))
            eod.append(round(data[trial_number][5]))
            electrode_frequency.append(round(data[trial_number][6]))
            start_time.append(round(data[trial_number][8]))
            choice_time.append(round(data[trial_number][9]))

        dicitonaries.append((dictionary))
        datas.append(data)
        trial_numbers.append(trial_number)


    return dicitonaries, keys, datas, rewarded_electrode, chosen_electrode, amplitude_1, amplitude_2, eod, electrode_frequency, start_time, choice_time


def get_data_versuchsprotokoll(fish):
    temperatures = []
    eods = []
    dates = []

    for f in fish:
        parts = f.strip().split(',')  #f.split(',') divides the data that is currently organised in lines at the comma
        # f.strip takes away the /n on every end of a line
        temperatures.append(float(parts[6]))
        eods.append(float(parts[4]))
        dates.append(parts[0])
    return temperatures, eods, dates


def get_file_names(fish):
    '''
    function makes a list of relacsfilenames and one list of videofilenames that are given in the 'versuchsprotokoll4.csv' file
    :param fish: the variable fish includes the lines of the 'versuchsprotokoll4.csv' file that are already ordered by fish
    :return: videofiles, relacs_files
    '''

    videofiles = []
    relacs_files = []
    fitting_dates = []


    for f in fish:
        parts = f.strip().split(',')  #f.split(',') divides the data that is currently organised in lines at the comma
        # f.strip takes away the /n on every end of a line
        if 'j' in parts[11] and 'j' in parts[13]: #at the moment only data is used, when relacsfile AND videofile is suitable
            videofiles.append(parts[10])
            relacs_files.append(parts[12])
            fitting_dates.append(parts[0])


    return videofiles, relacs_files, fitting_dates


def read_data(file):
    '''
    funktion liest 'versuchsprotokolle4.csv' datei und gibt die daten nach fischen sortiert zurueck
    :param file: ist der filename  (pfad) von 'versuchsprotokoll4.csv', in dem alle versuchsdurchlaeufe von versuch4 von jedem Fisch protokolliert sind
    :return: gibt einzelne lines der datei nach fischen sortiert zurueck
    '''

    lines = file.readlines()  # reads csv file and gives it back as lines, the columns of the file are divided through comma
    file.close()  # file hast to be closed after reading it
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

#code starts here
file = file('/home/plappert/Documents/Versuchsprotokoll_Versuch4.csv')
print file
chip, chap, alfons, trixi, krummschwanz, hermes = read_data(file)
videofiles1, relacs_files1, dates1 = get_file_names(chip)
videofiles2, relacs_files2, dates2 = get_file_names(chap)
videofiles3, relacs_files3, dates3 = get_file_names(alfons)
videofiles4, relacs_files4, dates4 = get_file_names(trixi)
videofiles5, relacs_files5, dates5 = get_file_names(krummschwanz)
videofiles6, relacs_files6, dates6 = get_file_names(hermes)


fish1 = '2015albi02'
fish2 = '2015albi01'
fish3 = '2014albi08'
fish4 = '2013albi14'
fish5 = '2013albi09'
fish6 = '2012albi01'

#temperatures1, eods1, dates1 = get_data_versuchsprotokoll(chip)
#temperatures4_2, eods4_2, dates4_2 = get_data_versuchsprotokoll(chap)
#temperatures4_3, eods4_3, dates4_3 = get_data_versuchsprotokoll(alfons)
#temperatures4_4, eods4_4, dates4_4 = get_data_versuchsprotokoll(trixi)
#temperatures4_5, eods4_5, dates4_5 = get_data_versuchsprotokoll(krummschwanz)
#temperatures4_6, eods4_6, dates4_6 = get_data_versuchsprotokoll(hermes)

dicitonaries1, keys1, data1, rewarded_electrode1, chosen_electrode1, amplitude_a1, amplitude_b1, eod1, electrode_frequency1, start_time1, choice_time1 = read_relacs_files(relacs_files1)
dicitonaries2, keys2, data2, rewarded_electrode2, chosen_electrode2, amplitude_a2, amplitude_b2, eod2, electrode_frequency2, start_time2, choice_time2 = read_relacs_files(relacs_files2)
dicitonaries3, keys3, data3, rewarded_electrode3, chosen_electrode3, amplitude_a3, amplitude_b3, eod3, electrode_frequency3, start_time3, choice_time3 = read_relacs_files(relacs_files3)
dicitonaries4, keys4, data4, rewarded_electrode4, chosen_electrode4, amplitude_a4, amplitude_b4, eod4, electrode_frequency4, start_time4, choice_time4 = read_relacs_files(relacs_files4)
dicitonaries5, keys5, data5, rewarded_electrode5, chosen_electrode5, amplitude_a5, amplitude_b5, eod5, electrode_frequency5, start_time5, choice_time5 = read_relacs_files(relacs_files5)
dicitonaries6, keys6, data6, rewarded_electrode6, chosen_electrode6, amplitude_a6, amplitude_b6, eod6, electrode_frequency6, start_time6, choice_time6 = read_relacs_files(relacs_files6)


percentage_right_choices1, percentage_wrong_choices1, percentage_no_choices1, percentage_fitting_dates1, right_choice1, wrong_choice1, no_choice1 = performance_bar_plot(rewarded_electrode1, chosen_electrode1, dates1, fish1)
percentage_right_choices2, percentage_wrong_choices2, percentage_no_choices2, percentage_fitting_dates2, right_choice2, wrong_choice2, no_choice2 = performance_bar_plot(rewarded_electrode2, chosen_electrode2, dates2, fish2)
#percentage_right_choices3, percentage_wrong_choices3, percentage_no_choices3, percentage_fitting_dates3, right_choice3, wrong_choice3, no_choice3 = performance_bar_plot(rewarded_electrode3, chosen_electrode3, dates_relacsfiles3, fish3)
percentage_right_choices4, percentage_wrong_choices4, percentage_no_choices4, percentage_fitting_dates4, right_choice4, wrong_choice4, no_choice4 = performance_bar_plot(rewarded_electrode4, chosen_electrode4, dates4, fish4)
percentage_right_choices5, percentage_wrong_choices5, percentage_no_choices5, percentage_fitting_dates5, right_choice5, wrong_choice5, no_choice5 = performance_bar_plot(rewarded_electrode5, chosen_electrode5, dates5, fish5)
percentage_right_choices6, percentage_wrong_choices6, percentage_no_choices6, percentage_fitting_dates6, right_choice6, wrong_choice6, no_choice6 = performance_bar_plot(rewarded_electrode6, chosen_electrode6, dates6, fish6)

#performance_line_plot(percentage_right_choices1, percentage_fitting_dates1, fish1)
#performance_line_plot(percentage_right_choices2, percentage_fitting_dates2, fish2)
#performance_line_plot(percentage_right_choices3, percentage_fitting_dates3, fish3)
#performance_line_plot(percentage_right_choices4, percentage_fitting_dates4, fish4)
#performance_line_plot(percentage_right_choices5, percentage_fitting_dates5, fish5)
#performance_line_plot(percentage_right_choices6, percentage_fitting_dates6, fish6)


mean_performance1 = performance_without_no_choice(right_choice1, wrong_choice1, no_choice1, dates1, fish1)
mean_performance2 = performance_without_no_choice(right_choice2, wrong_choice2, no_choice2, dates2, fish2)
#mean_performance3 = performance_without_no_choice(right_choice3, wrong_choice3, no_choice3, dates3, fish3)
mean_performance4 = performance_without_no_choice(right_choice4, wrong_choice4, no_choice4, dates4, fish4)
mean_performance5 = performance_without_no_choice(right_choice5, wrong_choice5, no_choice5, dates5, fish5)
mean_performance6 = performance_without_no_choice(right_choice6, wrong_choice6, no_choice6, dates6, fish6)

p_value1 = binomialtest(right_choice1, wrong_choice1, no_choice1)
p_value2 = binomialtest(right_choice2, wrong_choice2, no_choice2)
#p_value3 = binomialtest(right_choice3, wrong_choice3, no_choice3)
p_value4 = binomialtest(right_choice4, wrong_choice4, no_choice4)
p_value5 = binomialtest(right_choice5, wrong_choice5, no_choice5)
p_value6 = binomialtest(right_choice6, wrong_choice6, no_choice6)