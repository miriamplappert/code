__author__ = 'plappert'

import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from IPython import embed


def eichhornia_bar_plot(eichhornia, eichhornia2, species):
    species_eichhornia = OrderedDict()
    for s, t in zip(species,eichhornia):  #times and dates get ordered in a dictionary. therefore the date is the key on that the times of the date can be accesed
        if s not in species_eichhornia.keys():
            species_eichhornia[s] = []
        species_eichhornia[s].append(t)
    eichhornia_eigenmannia_positiv = []
    eichhornia_eigenmannia_negativ = []
    eichhornia_pulsfish_positiv = []
    eichhornia_pulsfish_negativ = []
    eichhornia_nofish_positiv = []
    eichhornia_nofish_negativ = []


    for h in species_eichhornia['E']:
        if h is 'j':
            eichhornia_eigenmannia_positiv.append(h)
        if h is 'n':
            eichhornia_eigenmannia_negativ.append(h)

    for h in species_eichhornia['P']:
        if h is 'j':
            eichhornia_pulsfish_positiv.append(h)
        if h is 'n':
            eichhornia_pulsfish_negativ.append(h)

    for e in eichhornia2:
        if e is 'j':
            eichhornia_nofish_positiv.append(e)
        if e is 'n':
            eichhornia_nofish_negativ.append(e)

    grundwert_eigenmannia = float(len(eichhornia_eigenmannia_positiv)) + float(len(eichhornia_eigenmannia_negativ))
    grundwert_pulsfish = float(len(eichhornia_pulsfish_positiv)) + float(len(eichhornia_pulsfish_negativ))
    grundwert_nofish = float(len(eichhornia2))

    p_eigenmannia_positiv = float(len(eichhornia_eigenmannia_positiv)) / grundwert_eigenmannia * 100
    p_eigenmannia_negativ = float(len(eichhornia_eigenmannia_negativ)) / grundwert_eigenmannia * 100
    p_pulsfish_positiv = float(len(eichhornia_pulsfish_positiv)) / grundwert_pulsfish * 100
    p_pulsfish_negativ = float(len(eichhornia_pulsfish_negativ)) / grundwert_pulsfish * 100
    p_nofish_positiv = float(len(eichhornia_nofish_positiv)) / grundwert_nofish * 100
    p_nofish_negativ = float(len(eichhornia_nofish_negativ)) / grundwert_nofish * 100

    print len(eichhornia_eigenmannia_positiv)
    print grundwert_eigenmannia
    print p_nofish_negativ
    print p_eigenmannia_positiv
    print p_pulsfish_positiv


    N = 1
    ind = np.arange(N)  # the x locations for the groups
    width = 0.01  # the width of the bars
    fig, ax = plt.subplots()

    ## the bars
    first_bar = ax.bar(ind, p_eigenmannia_positiv, width, color='green')  # creates successful bar in green
    second_bar = ax.bar(ind + width, p_eigenmannia_negativ, width, color='red')  # creates unsuccessful bar in red
    third_bar = ax.bar(ind + 3*width, p_pulsfish_positiv, width, color='green')
    fourth_bar = ax.bar(ind + 4*width, p_pulsfish_negativ, width, color='red')
    fifth_bar = ax.bar(ind + 6*width, p_nofish_positiv, width, color='green')
    sixth_bar = ax.bar(ind + 7*width, p_nofish_negativ, width, color='red')

    ax.set_ylabel('Auftreten von Eichhornia [%]')
    ax.set_title('Eichhornia')
    xTickMarks = ['Eigenmannia', 'Pulsfisch', 'kein Fund']
    ax.set_xticks([ind+width, ind + 4*width, ind+7*width])
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=10)
    plt.ylim(0, 100)
    ax.legend((first_bar[0], second_bar[0]), ('Eichhornia vorhanden', 'Eichhornia nicht vorhanden'))
    plt.subplots_adjust(left=0.12, bottom=0.15)
    plt.savefig('Eichhornia_plot.pdf')
    plt.show()




def plants_bar_plot(plants, plants2, species):
    species_plants = OrderedDict()
    for s, t in zip(species,plants):  #times and dates get ordered in a dictionary. therefore the date is the key on that the times of the date can be accesed
        if s not in species_plants.keys():
            species_plants[s] = []
        species_plants[s].append(t)
    mean_plants = []
    std_plants= []
    for k in species_plants.keys():
        mean_plants.append(np.mean(species_plants[k]))
        std_plants.append(np.std(species_plants[k]))  #np.std is the command for standard deviation

    mean_plants_nofish = np.mean(plants2)
    std_plants_nofish = np.std(plants2)


    N = 1
    ind = np.arange(N)  # the x locations for the groups
    width = 0.01  # the width of the bars
    fig, ax = plt.subplots()

    ## the bars
    first_bar = ax.bar(ind, mean_plants[0], width, color='green', yerr=std_plants[0], error_kw=dict(elinewidth=2,ecolor='black'))  # creates successful bar in green
    second_bar = ax.bar(ind + 1.5*width, mean_plants[1], width, color='yellow',  yerr=std_plants[1], error_kw=dict(elinewidth=2,ecolor='black'))  # creates unsuccessful bar in red
    third_bar = ax.bar(ind + 3*width, mean_plants_nofish, width, color='red',  yerr=std_plants_nofish, error_kw=dict(elinewidth=2,ecolor='black'))
    ax.set_ylabel('Grad des Bewuchses (1-4)')
    ax.set_title('Bewuchs')
    ax.set_xticks(ind + width)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.ylim(0, 4.5)
    ax.legend((first_bar[0], second_bar[0],third_bar[0]), ('Eigenmannia', 'Pulsfisch', 'kein Fund'))
    plt.savefig('Bewuchs_plot.pdf')
    plt.show()


def conductivity_bar_plot(conductivity, conductivity2, species):
    species_conductivity = OrderedDict()
    for s, t in zip(species,conductivity):  #times and dates get ordered in a dictionary. therefore the date is the key on that the times of the date can be accesed
        if s not in species_conductivity.keys():
            species_conductivity[s] = []
        species_conductivity[s].append(t)
    mean_conductivities = []
    std_conductivities = []
    for k in species_conductivity.keys():
        mean_conductivities.append(np.mean(species_conductivity[k]))
        std_conductivities.append(np.std(species_conductivity[k]))  #np.std is the command for standard deviation

    mean_conductivity_nofish = np.mean(conductivity2)
    std_conductivity_nofish = np.std(conductivity2)
    

    N = 1
    ind = np.arange(N)  # the x locations for the groups
    width = 0.01  # the width of the bars
    fig, ax = plt.subplots()

    ## the bars
    first_bar = ax.bar(ind, mean_conductivities[0], width, color='green', yerr=std_conductivities[0], error_kw=dict(elinewidth=2,ecolor='black'))  # creates successful bar in green
    second_bar = ax.bar(ind + 1.5*width, mean_conductivities[1], width, color='yellow',  yerr=std_conductivities[1], error_kw=dict(elinewidth=2,ecolor='black'))  # creates unsuccessful bar in red
    third_bar = ax.bar(ind + 3*width, mean_conductivity_nofish, width, color='red',  yerr=std_conductivity_nofish, error_kw=dict(elinewidth=2,ecolor='black'))
    ax.set_ylabel('Leitfaehigkeit [Microsiemens pro cm]')
    ax.set_title('Leitfaehigkeit')
    ax.set_xticks(ind + width)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.ylim(0, 130)
    ax.legend((first_bar[0], second_bar[0],third_bar[0]), ('Eigenmannia', 'Pulsfisch', 'kein Fund'))
    plt.savefig('Leitfaehigkeit_plot.pdf')
    plt.show()
    


def temperature_bar_plot(temperature, temperature2, species):
    species_temperature = OrderedDict()
    for s, t in zip(species,temperature):  #times and dates get ordered in a dictionary. therefore the date is the key on that the times of the date can be accesed
        if s not in species_temperature.keys():
            species_temperature[s] = []
        species_temperature[s].append(t)
    mean_temperatures = []
    std_temperatures = []
    for k in species_temperature.keys():
        mean_temperatures.append(np.mean(species_temperature[k]))
        std_temperatures.append(np.std(species_temperature[k]))  #np.std is the command for standard deviation

    mean_temperature_nofish = np.mean(temperature2)
    std_temperature_nofish = np.std(temperature2)
    print species
    print temperature
    print temperature2
    print species_temperature
    print std_temperatures

    N = 1
    ind = np.arange(N)  # the x locations for the groups
    width = 0.01  # the width of the bars
    fig, ax = plt.subplots()

    ## the bars
    first_bar = ax.bar(ind, mean_temperatures[0], width, color='green', yerr=std_temperatures[0], error_kw=dict(elinewidth=2,ecolor='black'))  # creates successful bar in green
    second_bar = ax.bar(ind + 1.5*width, mean_temperatures[1], width, color='yellow',  yerr=std_temperatures[1], error_kw=dict(elinewidth=2,ecolor='black'))  # creates unsuccessful bar in red
    third_bar = ax.bar(ind + 3*width, mean_temperature_nofish, width, color='red',  yerr=std_temperature_nofish, error_kw=dict(elinewidth=2,ecolor='black'))
    ax.set_ylabel('Temperatur [Grad Celsius]')
    ax.set_title('Temperatur')
    ax.set_xticks(ind + width)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.ylim(0, 40)
    ax.legend((first_bar[0], second_bar[0],third_bar[0]), ('Eigenmannia', 'Pulsfisch', 'kein Fund'))
    plt.savefig('Temperatur_plot.pdf')
    plt.show()



def read_data(file):

    lines = file.readlines()  # reads csv file and gives it back as lines, the colunms of the file are divided through comma
    file.close()  # file hast to be closed after reading it

    fish=[]
    no_fish=[]

    for l in lines:
        if 'ja' in l:
            fish.append(l)
        if 'nein' in l:
            no_fish.append(l)
    return fish, no_fish

def order_data(list):
    plants = []  # empty list named after the first fish
    conductivity= []
    temperature = []
    species = []
    water_size = []
    floating = []
    eichhornia = []

    for f in list:
        parts = f.strip().split(',')  #f.split(',') devides the data, that is currently organised in lines, at the comma
        # f.strip takes away the /n on every end of a line
        plants.append(float(parts[0]))  # command appends the first element of the line to the list date
        if float(parts[1]) == 0:
            continue
        else:
            conductivity.append(float(parts[1]))  #float converts strings to numbers
        if float(parts[2]) == 0:
            continue
        else:
            temperature.append(float(parts[2]))
        species.append(parts[3])
        water_size.append(float(parts[4]))
        floating.append(parts[5])
        eichhornia.append(parts[7])

    return plants, conductivity, temperature, species, water_size, floating, eichhornia




if __name__ == '__main__':  # the code doesn't run if someone is adding it to his or her code it only runs if it is the main code
    file = file('/home/plappert/Documents/brasilien_auswertung.csv')  #the file that we need gets the variable file
    print file


fish, no_fish = read_data(file)
plants, conductivity, temperature, species, water_size, floating, eichhornia = order_data(fish)
plants2, conductivity2, temperature2, species2, water_size2, floating2, eichhornia2 = order_data(no_fish)

temperature_bar_plot(temperature,temperature2, species)
conductivity_bar_plot(conductivity,conductivity2,species)
plants_bar_plot(plants, plants2, species)
eichhornia_bar_plot(eichhornia, eichhornia2, species)