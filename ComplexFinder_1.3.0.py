# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 13:36:15 2021

@author: Marios
"""
"""
ComplexFinder_1.3.0

The current version of the algorithm does the following:
    
    1) Reads the PPI datafiles and returns an appropriate data structure.
    2) Reads the HGNC datafiles and returns an appropriate datastructure.
    3) Compares them and visualizes the results as bar charts.
    
"""

def ppi_reader(datafile):
   """
   1) Reads the tsv PPI input files.
   2) Stores each line in a list object.
   
   Returns
   -------
   datatable (list object)

   """
   interactions = 0 # variable for interactions in each datafile
       
   # Open data file in "read" mode.
   filehandle = open(datafile, 'r')
   headers = filehandle.readline()
   
   datatable = []
   
   # The interactions (lines) in the data file are calculated.
   # Lines are splitted and data are stored in a data structure called datatable.
   for line in filehandle:
       interactions += 1
       line.strip()
       data = line.split('\t')
       datatable.append(data)
   
   return(datatable)

def complex_reader(datafile):
    """
    1) Reads the HGNC input file.
    2) Creates a dictionary (key = complex name, value = list of proteins).
    
    Parameters
    ----------
    datafile: list object returned by ppi_reader() function. 
       
    Returns
    -------
    complex_dict (dictionary object)

    """
       
    # Open data file in "read" mode.
    filehandle = open(datafile, 'r')
    headers = filehandle.readline()
    
    datatable = []
    complex_dict = {}
    
    # The interactions (lines) in the data file are calculated.
    # Lines are splitted and data are stored in a data structure called datatable.
    for line in filehandle:
        line.strip()
        data = line.split('\t')
        datatable.append(data)
        
    for entry in datatable:
        if entry[0] not in complex_dict.keys():
            complex_dict.update({entry[0]:[entry[2]]})
        else:
            complex_dict[entry[0]].append(entry[2])
   
    return(complex_dict)

def complex_finder(oct4_datatable, e7_datatable, complex_dictionary):
    """
    1) Reads the PPI and HGNC data structures and compares them.
    2) Calls the plotter() fuction to visualize the results as bar charts.

    Parameters
    ----------
    oct4_datatable: list object returned by ppi_reader() function.
    e7_datatable: list object returned by ppi_reader() function.
    complex_dictionary: dictionary returned by complex_reader() function.
    
    Returns
    -------
    Bar Charts

    """
    
    oct4_ppis = [] # list object where Oct4 protein interactors are saved
    e7_ppis = [] # list object where E7 protein interactors are saved
    
    for entry in oct4_datatable:
        oct4_ppis.append(entry[5])
    
    for entry in e7_datatable:
        e7_ppis.append(entry[5])
    
    oct4_dict = {} # dictionary (key = protein complex, value = oct4_list)
    e7_dict = {} # dictionary (key = protein complex, value = e7_list)
    both_dict = {} # dictionary (key = protein complex, value = both_list)
    
    for complexx, component in complex_dictionary.items():
        # lists where protein interactors found in protein complexes are saved
        oct4_list = [] # Oct4 protein interactors found in complexes
        e7_list = [] # E7 protein interactors found in complexes
        both_list = [] # Common protein interactors found in complexes
        
        for prot in component:
            if prot in oct4_ppis:
                oct4_list.append(prot)
                oct4_dict.update({complexx:oct4_list})
            if prot in e7_ppis:
                e7_list.append(prot)
                e7_dict.update({complexx:e7_list})
            if prot in e7_ppis:
                if prot in oct4_ppis:
                    both_list.append(prot)
                    both_dict.update({complexx:both_list})
    
    # Calls to plotter() function for the visualization of the results.
    plotter(oct4_dict, complex_dictionary, 'Complexes with proteins interacting with Oct4.', 8)
    plotter(e7_dict, complex_dictionary, 'Complexes with proteins interacting with E7.', 10)
    plotter(both_dict, complex_dictionary, 'Complexes with proteins interacting with Oct4 & E7.', 12)
    exporter(oct4_dict, 'Supplementary Table 8.txt')
    exporter(e7_dict, 'Supplementary Table 9.txt')
    exporter(both_dict, 'Supplementary Table 10.txt')

def plotter(dictionary, complex_dictionary, title, font):
    """
    1) Utilizes matplotlib library to visualize the results of the comparative
    analysis.
    *** plotter() function is called inside the complex_finder() function.

    Parameters
    ----------
    dictionary: dictionary created in complex_finder() function.
    complex_dictionary: dictionary returned by complex_reader() function.
    title: title of bar chart
    font: numeric value of letter fontsize.

    Returns
    -------
    Bar Charts.

    """
    # import of necessary modules
    import matplotlib.pyplot as plt
    import numpy as np
    
    percentages = {} # key = complex name, value = interactors percentage in complex
    fractions = {} # key = complex name, value = No of proteins in complex & No
                # of interactors in complex
    
    for complexx, component in dictionary.items():
        percentage = (len(component)/len(complex_dictionary[complexx])*100)
        r_percentage = round(percentage,3)
        percentages.update({complexx:r_percentage})
        fractions.update({complexx:[str(len(component)), str(len(complex_dictionary[complexx]))]})
        
    # plot settings
    plt.rcdefaults()
    fig, ax = plt.subplots(figsize=(10,7))
    
    # per dictionary is sorted decreasingly based on the percentages of protein
    # interactors in complexes and the names of the sorted complexes are saved
    # in the sorted labels list.
    s_percentages = sorted(percentages.items(),key = lambda x:x[1], reverse=True)
    s_fractions = []
    
    for entry in s_percentages:
        s_fractions.append(fractions[entry[0]])
   
    
    complex_list = [] # list where the sorted complexes are saved
    for i in s_percentages:
        complex_list.append(i[0])
    
    y_pos = np.arange(len(complex_list)) # position of each complex on the plot
    
    percentage = [] # list where the sorted percentages are saved
    for i in s_percentages:
       percentage.append(i[1])

    # list where the sorted percentages and fractions are saved
    percentage_fraction = []
    for p,f in zip(percentage,s_fractions):
       joined_f = '/'.join(f)
       percentage_fraction.append([p,joined_f])
    
    # Remove axes splines
    for s in ['top','bottom','left','right']:
        ax.spines[s].set_visible(False)
    
    # Remove x,y Ticks
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    
    # Add padding between axes and labels
    ax.xaxis.set_tick_params(pad=5)
    ax.yaxis.set_tick_params(pad=10)
    
    # Add x,y gridlines
    ax.grid(b=True, color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)

    ax.barh(y_pos, percentage, align='center')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(complex_list)
    ax.invert_yaxis()  # labels read top-to-bottom
    
    for i, v in enumerate(percentage_fraction):
        ax.text(x=v[0] + 1, y=i, s=v[1], color='black', fontsize=font, ha='left', va='center')
        
    ax.set_xlabel('Percentage (%)') # x-axis label
    ax.set_title(title) # plot title
    
    plt.show()

def exporter(dictionary, out_file):
    """
    1) Takes as argument the dictionary object cretaed in complex_finder().
    2) Outputs a tsv file with two columns (Complex, Proteins).

    Parameters
    ----------
    dictionary : dictionary object created in complex_finder()

    Returns
    -------
    tsv output file

    """
    import csv
    
   # A tsv data file with interactors is created.
    with open(out_file,'w', newline = '') as output_file:
        tsvfile = csv.writer(output_file, delimiter = '\t')
        tsvfile.writerow(['Complex', 'Proteins'])  
        
        for key, value in dictionary.items():
            tsvfile.writerow([key, value])

import os
import sys

sys.path.append(os.path.realpath('..'))
dirname = os.path.dirname(__file__)

all_oct4 = os.path.join(dirname, '..\Supplementary File 2\PPI data\All_Oct4.txt')
all_e7 = os.path.join(dirname, '..\Supplementary File 2\PPI data\All_E7.txt')
hs_oct4 = os.path.join(dirname, '..\Supplementary File 2\PPI data\Oct4_Homo_sapiens.txt')
hpv16_e7 = os.path.join(dirname, '..\Supplementary File 2\PPI data\E7_HPV16.txt')
hgnc = os.path.join(dirname, '..\Supplementary File 2\HGNC data\complexes.txt')

oct4_datatable = ppi_reader(all_oct4)
e7_datatable = ppi_reader(all_e7)
complexes = complex_reader(hgnc)
complex_finder(oct4_datatable,e7_datatable,complexes)