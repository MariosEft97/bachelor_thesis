import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import ast

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

def chip_data_reader(directory: str) -> list:
   
    """
    The function reads the tsv PPI ChIP-Atlas input files and stores them in a pandas DataFrame.
    
    Parameters:
        directory (str): path to data directory
    
    Returns: 
        df_list (list): data structure to save data
    """
    
    if not isinstance(directory, str):
        raise TypeError("directory must be specified as a string")
    
    else:
        
        df_list = []
        
        # iterate over files in
        # that directory
        for filename in os.listdir(directory):
            f = os.path.join(directory, filename)
            # checking if it is a file
            if os.path.isfile(f):
                df = pd.read_csv(f, sep = '\t', header = 0)
                df_list.append(df)

        return(df_list)

# ----------------------------------------------------------------------------------------------------

def deg_reader(directory: str) -> list:
   
    """
    The function reads the TCGA input files and stores each line as a list in a list object.
    
    Parameters:
        directory (str): path to data directory
    
    Returns: 
        deg_data_list (list): data structure to save data
    """
    
    if not isinstance(directory, str):
        raise TypeError("directory must be specified as a string")
    
    else:
    
        degs = 0 # variable for differentially expressed genes in each datafile
        deg_data_list = []
        
        # iterate over files in
        # that directory
        for filename in os.listdir(directory):
            f = os.path.join(directory, filename)
            # checking if it is a file
            if os.path.isfile(f):

                # Open data file in "read" mode.
                filehandle = open(f, 'r')

                # Input files are in dictionary format.
                # ast library is used for reading such files.
                contents = filehandle.read()
                tcga_dict = ast.literal_eval(contents)
    
                filehandle.close()
    
                deg = tcga_dict['associations']

                # List object for storing gene name and expression status.
                # 1 = High expression
                # -1 = Low expression
                deg_data = []
    
                for line in deg:
                    gene = line['gene']['symbol']
                    expression = line['thresholdValue']
                    deg_data.append([gene, expression])

                # Number of DEGs is calculated.
                for line in deg_data:
                    degs+=1
                
                deg_data_list.append(deg_data)
            
        return(deg_data_list)

# ----------------------------------------------------------------------------------------------------

def chip_df_constructor(input_df_list: list, protein: str) -> pd.DataFrame:
    """
    The function takes as argument the data structure returned from the chip_data_reader() function.
    It filters the gene targets and return a DataFrame with 4 columns (TF, Target Gene, macs2, q-value).

    Parameters:
        input_df_list (list): data structure returned from the chip_data_reader() function
        protein (str): POU5F1 or E7
    
    Returns:
        df (pd.DataFrame): data structure to save results

    """
    
    if not isinstance(input_df_list, list):
        raise TypeError("input_df_list must be specified as a list")
    
    elif not isinstance(protein, str):
        raise TypeError("protein must be specified as a string\noptions: POU5F1 or E7")
    
    else:
        
        output_df_list = []
        
        for file in input_df_list:
            file.drop(file.columns[[1,-1]],inplace=True, axis=1)

            headers = list(file.columns.values)
            headers.pop(0)

            header_dict = {}

            for header in headers:
                header_cell = header.split('|')
                experiment = header_cell[0]
                cell_line =header_cell[1]

                if cell_line not in header_dict.keys():
                    header_dict.update({cell_line:[experiment]})
                else:
                    header_dict[cell_line].append(experiment)

            inside_loop_df_list = []

            for cell_line, experiments in header_dict.items():

                header_list = []

                for experiment in experiments:
                    header = experiment + '|' + cell_line
                    header_list.append(header)

                header_list.insert(0,'Target_genes')

                df = file.loc[:, header_list]
                df['Average'] = df.mean(axis = 1)
                df_0 = df[df['Average']>0]
                inside_loop_df_list.append(df_0)


            target_average_dict = {}

            for df in inside_loop_df_list:

                new_df = df.loc[:, ['Target_genes', 'Average']]

                target_genes = new_df['Target_genes'].tolist()

                average = new_df['Average'].tolist()

                for gene, avrg in zip(target_genes, average):

                    occurance = 1

                    if gene not in target_average_dict.keys():

                        target_average_dict.update({gene:[avrg, occurance]})

                    else:

                        target_average_dict[gene][0] += avrg
                        target_average_dict[gene][1] += 1

            final_target_average_dict = {}

            for gene in target_average_dict.keys():

                macs2 = target_average_dict[gene][0]*target_average_dict[gene][1]/len(inside_loop_df_list)
                final_target_average_dict.update({gene:macs2})


            final_df = pd.DataFrame(final_target_average_dict.items(), columns = ['Target_genes', 'Average'])
            new_avg = final_df['Average'].tolist()
            tar = final_df['Target_genes'].tolist()

            test_list = []

            for (x, y) in zip(tar, new_avg):
                test_list.append((x, y))

            n = int(5*len(test_list)/100)

            res = sorted(test_list, key = lambda x: x[1], reverse = True)[:n]

            target = []
            macs2 = []
            for i in res:
                target.append(i[0])
                macs2.append(i[1])

            qval = []

            for score in macs2:
                val = float(10**(-score/10))
                qval.append(val)

            up_target = []
            for tar in target:
                up_target.append(tar.upper())

            source = [protein]*len(up_target) # TF column

            df = pd.DataFrame({'TF': source,
                                'Target Gene': up_target,
                                'macs2': macs2,
                                'q-value': qval})
            df.index += 1
            
            output_df_list.append(df)
        
        all_df = pd.concat(output_df_list, axis=0)
        
        datatable = [] # list object where the dataframe data are stored
        tgdict = {} # key = gene name, value = [mac2, q-value, occurance]
        ftgdict = {} # key = gene name, value = [overall macs2, overall q-value]
    
        all_df.drop('TF', inplace=True, axis=1)
    
        for i in range(len(all_df)):
            datatable.append([all_df.iloc[i, 0], all_df.iloc[i, 1], all_df.iloc[i, 2]])

        for entry in datatable:
            gene = entry[0]
            macs2 = entry[1]
            qval = entry[2]
            occurance = 1

            if gene not in tgdict.keys():
                tgdict.update({gene:[macs2, qval, occurance]})
            else:
                tgdict[gene][2] += 1
                tgdict[gene][0] += macs2
                tgdict[gene][1] += qval

        for gene, val in tgdict.items():
            fmacs = tgdict[gene][0]/tgdict[gene][2]
            fqval = tgdict[gene][1]/tgdict[gene][2]
            ftgdict.update({gene:[fmacs, fqval]})


        gene = ftgdict.keys() # Target gene column
        source = [protein]*len(gene) # TF column
        values = ftgdict.values()
        macs2 = [] # mac2 column
        qval = [] # q-value column

        for i in values:
            macs2.append(i[0])
            qval.append(i[1])

        df = pd.DataFrame({'TF': source,
                           'Target Gene': gene,
                           'macs2': macs2,
                           'q-value': qval})

        df.index += 1        
        return(df)
    
# ----------------------------------------------------------------------------------------------------

def deg_df_constructor(deg_data_list: list) -> pd.DataFrame:
    """
    The function takes as an argument the data structure returned from deg_reader() function.
    It reads gene targets and returns a DataFrame with 2 columns (DEG, Expression).

    Parameters:
        deg_data_list (list): data structure returned from the deg_reader() function
    
    Returns:
        df (pd.DataFrame): data structure to save results

    """
    
    if not isinstance(deg_data_list, list):
        raise TypeError("deg_data_list must be specified as a list")
    
    else:
    
        output_df_list = []

        for datatable in deg_data_list:

            genes = [] # list where DEGs are stored
            high_low = [] # list where expression values are stored (-1,1)

            for entry in datatable:
                genes.append(entry[0])
                high_low.append(entry[1])

            df = pd.DataFrame({'Genes': genes,
                               'Expression': high_low})

            output_df_list.append(df)

        all_df = pd.concat(output_df_list)
        
        datatable = [] # list object where the dataframe data are stored 
        lowdict = {} # key = low expressed gene name, value = (expression, occurance)
        highdict = {} # key = high expressed gene name, value = (expression, occurance)
        flowdegdict = {} # filtered dictionary based on gene occurance
        fhighdegdict = {} # filtered dictionary based on gene occurance
    
        for i in range(len(all_df)):
            datatable.append([all_df.iloc[i, 0], all_df.iloc[i, 1]])

        for entry in datatable:
            gene = entry[0]
            expression = entry[1]
            occurance = 1
        
            if expression == -1:
                if gene not in lowdict.keys():
                    lowdict.update({gene:[expression, occurance]})
                else:
                    lowdict[gene][1] += 1
            if expression == 1:
                if gene not in highdict.keys():
                    highdict.update({gene:[expression, occurance]})
                else:
                    highdict[gene][1] += 1

        lowfiltered = []
        highfiltered = []

        # Low expressed genes are not filterd as their numbers are low compared to
        # high expressed genes.
        for gene, val in lowdict.items():
            lowfiltered.append(gene)
            flowdegdict.update({gene:val})

        # High expressed genes are filtered by their occurance.
        # Only genes found in all conditions are stored.
        for gene, val in highdict.items():
            if val[1]>1:
                highfiltered.append(gene)
                fhighdegdict.update({gene:val})

        lowgenes = flowdegdict.keys()
        lowexpression = flowdegdict.values()

        highgenes = fhighdegdict.keys()
        highexpression = fhighdegdict.values()

        genes = list(highgenes) + list(lowgenes)
        expression = list(highexpression) + list(lowexpression)

        df = pd.DataFrame({'DEGs': genes,
                            'Expression': expression})

        return(df)

# ----------------------------------------------------------------------------------------------------

def gene_targets(merged_chip_df: pd.DataFrame, merged_deg_df: pd.DataFrame, protein: str) -> pd.DataFrame:
    """
    The function takes as arguments the data structures returned by chip_df_constructor()
       and deg_df_constructor() functions, counts the occurance of genes and returns
       a merged DataFrame with 4 columns (TF, Target Gene, macs2, q-value).

    Parameters:
        merged_chip_df (pd.DataFrame): data structure returned from chip_df_constructor()
        merged_deg_df (pd.DataFrame): data structure returned from deg_df_constructor()
        protein (str): POU5F1 or E7 or POU5F1/E7

    Returns:
        df (pd.DataFrame): data structure to save results
    """
     
    if not isinstance(merged_chip_df, pd.DataFrame):
        raise TypeError("merged_chip_df must be specified as a pandas DataFrame")
    
    elif not isinstance(merged_deg_df, pd.DataFrame):
        raise TypeError("merged_deg_df must be specified as a pandas DataFrame")
    
    elif not isinstance(protein, str):
        raise TypeError("protein must be specified as a string\noptions: POU5F1 or E7 or POU5F1/E7")
    
    else:
    
        chip_targets = merged_chip_df['Target Gene'].tolist()
        q_value = merged_chip_df['q-value'].tolist()
        deg_genes = merged_deg_df['DEGs'].tolist()
        expression = merged_deg_df['Expression'].tolist()
    
        chip_dict = {chip_targets[i]: q_value[i] for i in range(len(chip_targets))}
        deg_dict = {deg_genes[i]: expression[i] for i in range(len(deg_genes))}
    
        commons = {}
        for gene in deg_dict.keys():
            if gene in chip_dict.keys():
                commons.update({gene:[chip_dict[gene],deg_dict[gene]]})
    
        genes = [] # list of TCGA genes that are probably targeted by Oct4
        qval = [] # list of the q-values
        exp = [] # list of the expression status (1 or -1)
        exp_1 = [] # list of the expression status (High or Low)
        for gene, value in commons.items():
            genes.append(gene)
            qval.append(value[0])
            exp.append(value[1][0])

        for i in exp:
            if i == 1.0:
                exp_1.append('High')
            else:
                exp_1.append('Low')

        source = [protein]*len(genes) 

        df = pd.DataFrame({'TF': source,
                           'Genes':genes,
                           'q-value':qval,
                           'Expression':exp_1})
                    
        return(df)

# ----------------------------------------------------------------------------------------------------

def common_gene_targets(gene_targets_a: pd.DataFrame, gene_targets_b: pd.DataFrame) -> pd.DataFrame:
    
    """
    The function takes as arguments the DataFrame objects returned from gene_targets() function.
    It compares the two input DataFrames and returns new DataFrame containing only the common interactors.

    Parameters:
        gene_targets_a (pd.DataFrame): data structure returned from gene_targets() function
        gene_targets_b (pd.DataFrame): data structure returned from gene_targets() function

    Returns:
        df (pd.DataFrame): data structure containing only common gene targets

    """
    
    if not isinstance(gene_targets_a, pd.DataFrame):
        raise TypeError("gene_targets_a must be specified as a pandas DataFrame")
    
    elif not isinstance(gene_targets_b, pd.DataFrame):
        raise TypeError("gene_targets_b must be specified as a pandas DataFrame")
    
    else:
    
        common = []
        expression = []
        regulation = []
    
        for gene, exp in zip(list(gene_targets_a["Genes"]), list(gene_targets_a["Expression"])):
            if gene in list(gene_targets_b["Genes"]):
                common.append(gene)
                expression.append(exp)
    
        for exp in expression:
            if exp == 'High':
                regulation.append('activation')
            else:
                regulation.append('repression')
    
        source = ['POU5F1/E7']*len(common)
    
        df = pd.DataFrame({'Source': source,
                           'Target': common,
                           'Regulation': regulation,
                           'Expression': expression})
    
        return(df)

# ----------------------------------------------------------------------------------------------------