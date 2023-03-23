import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------------------------------------

def ppi_reader(data_filepath: str) -> list:
   
    """
    The function reads the tsv PPI input files and stores each line as a list in a list object.
    
    Parameters:
        data_filepath (str): filepath to data
    
    Returns: 
        datatable (list): data structure to save data
    """
    
    if not isinstance(data_filepath, str):
        raise TypeError("data_filepath must be specified as a string")
    
    else:
    
        interactions = 0 # variable for interactions in each datafile

        # Open data file in "read" mode.
        filehandle = open(data_filepath, 'r')
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

# ----------------------------------------------------------------------------------------------------

def complex_reader(data_filepath: str) -> dict:
    
    """
    The function reads the HGNC input file and stores the data in a
    dictionary (key = complex name, value = list of proteins).
    
    Parameters:
        data_filepath (str): filepath to data
    
    Returns: 
        complex_dict (dict): data structure to save data
    """
    
    if not isinstance(data_filepath, str):
        raise TypeError("data_filepath must be specified as a string")
    
    else:
       
        # Open data file in "read" mode.
        filehandle = open(data_filepath, 'r')
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

# ----------------------------------------------------------------------------------------------------

def complex_finder(datatable_a: list, datatable_b: list, complex_dictionary: dict) -> dict:
    
    """
    The function reads the PPI and HGNC data structures and compares them and
    visualizes the results as bar charts.

    Parameters:
        datatable_a (list): data structure returned from ppi_reader() function.
        datatable_b (list): data structure returned from ppi_reader() function.
        complex_dictionary (dictionary): data structure returned from complex_reader() function.
    
    Returns:
        a_dict (dict): data structure to save results
        b_dict (dict): data structure to save results
        both_dict (dict): data structure to save results
    """
    
    if not isinstance(datatable_a, list):
        raise TypeError("datatable_a must be specified as a list")
    
    elif not isinstance(datatable_b, list):
        raise TypeError("datatable_b must be specified as a list")
        
    elif not isinstance(complex_dictionary, dict):
        raise TypeError("complex_dictionary must be specified as a dictionary")
    
    else:
    
        ppis_a = [] # list object where Oct4 protein interactors are saved
        ppis_b = [] # list object where E7 protein interactors are saved

        for entry in datatable_a:
            ppis_a.append(entry[5])

        for entry in datatable_b:
            ppis_b.append(entry[5])

        a_dict = {} # dictionary (key = protein complex, value = oct4_list)
        b_dict = {} # dictionary (key = protein complex, value = e7_list)
        both_dict = {} # dictionary (key = protein complex, value = both_list)

        for complexx, component in complex_dictionary.items():
            # lists where protein interactors found in protein complexes are saved
            oct4_list = [] # Oct4 protein interactors found in complexes
            e7_list = [] # E7 protein interactors found in complexes
            both_list = [] # Common protein interactors found in complexes

            for prot in component:
                if prot in ppis_a:
                    oct4_list.append(prot)
                    a_dict.update({complexx:oct4_list})
                if prot in ppis_b:
                    e7_list.append(prot)
                    b_dict.update({complexx:e7_list})
                if prot in ppis_b:
                    if prot in ppis_a:
                        both_list.append(prot)
                        both_dict.update({complexx:both_list})

        
        return(a_dict, b_dict, both_dict)
    
# ----------------------------------------------------------------------------------------------------

def complex_plotter(complex_finder_dictionary: dict, complex_reader_dictionary: dict, title: str, font: int) -> None:
    
    """
    The function uses matplotlib to visualize the results of complex_finder() function.

    Parameters:
        complex_finder_dictionary (dict): data structured returned from complex_finder() function
        complex_reader_dictionary (dict): data structured returned from complex_reader() function
        title (str): title of bar chart
        font (integer): numeric value of letter fontsize

    Returns:
        None
    """
    
    if not isinstance(complex_finder_dictionary, dict):
        raise TypeError("complex_finder_dictionary must be specified as a dictionary")
    
    elif not isinstance(complex_reader_dictionary, dict):
        raise TypeError("complex_reader_dictionary must be specified as a dictionary")
        
    elif not isinstance(title, str):
        raise TypeError("title must be specified as a string")
    
    elif not isinstance(font, int):
        raise TypeError("font must be specified as an integer value")
    
    else:
        
        percentages = {} # key = complex name, value = interactors percentage in complex
        fractions = {} # key = complex name, value = No of proteins in complex & No of interactors in complex
    
        for complexx, component in complex_finder_dictionary.items():
            percentage = (len(component)/len(complex_reader_dictionary[complexx])*100)
            r_percentage = round(percentage,3)
            percentages.update({complexx:r_percentage})
            fractions.update({complexx:[str(len(component)), str(len(complex_reader_dictionary[complexx]))]})
        
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
        ax.grid(color='grey', linestyle='-.', linewidth=0.5, alpha=0.2)

        ax.barh(y_pos, percentage, align='center')
        ax.set_yticks(y_pos)
        ax.set_yticklabels(complex_list)
        ax.invert_yaxis()  # labels read top-to-bottom

        for i, v in enumerate(percentage_fraction):
            ax.text(x=v[0] + 1, y=i, s=v[1], color='black', fontsize=font, ha='left', va='center')

        ax.set_xlabel('Percentage (%)') # x-axis label
        ax.set_title(title) # plot title

        plt.show()
        
        return(None)

# ----------------------------------------------------------------------------------------------------