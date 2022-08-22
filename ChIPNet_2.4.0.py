# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 17:22:11 2021

@author: Marios
"""

"""
ChIPNet_2.4.0.py

The current version of the algorithm does the following:
    
    1) Reads the PPI, ChIP-Atlas and TCGA input datafiles.
    2) Returns an appropriate data structure for each input file.
    3) Compares them and returns the predicted gene targets.
    4) Visualizes results in Cytoscape platform.
"""

def ppi_reader(datafile):
    """
    1) Reads the tsv PPI input files.
    2) Stores each line in a list object.
    
    Parameters
    ----------
    datafile (tsv PPI data files)

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
    
    print(f'PPIs = {interactions}.')
    return(datatable)

def oct4_chip_reader(datafile):
    """
    1) Reads the Oct4 ChIP-Atlas input files.
    2) Stores each line in a pandas Dataframe.
    
    Parameters
    ----------
    datafile (tsv ChIP-Atlas data files)

    Returns
    -------
    file (pandas Dataframe)

    """
    import pandas as pd
    
    file = pd.read_csv(datafile, sep = '\t', header = 0)
    #print(file.info)
    
    return(file)

def e7_chip_reader(datafile):
    """
    1) Reads the E7 ChIP-Atlas input files.
    2) Stores each line in a pandas Dataframe.
    
    Parameters
    ----------
    datafile (tsv ChIP-Atlas data files)

    Returns
    -------
    file (pandas Dataframe)

    """    
    import pandas as pd
    
    file = pd.read_csv(datafile, sep = '\t', header = 0)
    #print(file.info)
    return(file)

def deg_reader(datafile):
    """
    1) Reads the TCGA input files.
    2) Stores each line in a list object.

    Parameters
    ----------
    datafile (tsv TCGA data files)

    Returns
    -------
    deg_data (list object)

    """
    import ast
    
    degs = 0 # variable for differentially expressed genes in each datafile
        
    # Open data file in "read" mode.
    filehandle = open(datafile, 'r')
    
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
    
    print(f'DEGs = {degs}')
    return(deg_data)
 
def homo_oct4_chip_df_constructor(file):
    """
    1) Takes as argument list objects returned by the chip_reader() functions.
    2) Filters the gene targets.
    3) Returns DataFrames with 4 columns (TF, Target Gene, macs2, q-value).

    Parameters
    ----------
    chip_datatable (list object returned by chip_reader() functions)
    
    Returns
    -------
    pandas DataFrame

    """
    import pandas as pd
    
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
        
    df_list = []
    
    for cell_line, experiments in header_dict.items():
        
        header_list = []
        
        for experiment in experiments:
            header = experiment + '|' + cell_line
            header_list.append(header)
        
        header_list.insert(0,'Target_genes')
        
        df = file.loc[:, header_list]
        df['Average'] = df.mean(axis = 1)
        df_0 = df[df['Average']>0]
        df_list.append(df_0)
    
       
    target_average_dict = {}
    
    for df in df_list:
        
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
        
        macs2 = target_average_dict[gene][0]*target_average_dict[gene][1]/len(df_list)
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
    
    source = ['POU5F1']*len(up_target) # TF column
    
    df = pd.DataFrame({'TF': source,
                        'Target Gene': up_target,
                        'macs2': macs2,
                        'q-value': qval})
    df.index += 1
    print(f'Predicted Gene Targets = {len(df.index)}')
    #print(df)
    return(df)

def mus_oct4_chip_df_constructor(file):
    """
    1) Takes as argument list objects returned by the chip_reader() functions.
    2) Filters the gene targets.
    3) Returns DataFrames with 4 columns (TF, Target Gene, macs2, q-value).

    Parameters
    ----------
    chip_datatable (list object returned by chip_reader() functions)
    
    Returns
    -------
    pandas DataFrame

    """
    import pandas as pd
    
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
    
    df_list = []
    
    for cell_line, experiments in header_dict.items():
        
        header_list = []
        
        for experiment in experiments:
            header = experiment + '|' + cell_line
            header_list.append(header)
        
        header_list.insert(0,'Target_genes')
        
        df = file.loc[:, header_list]
        df['Average'] = df.mean(axis = 1)
        df_0 = df[df['Average']>0]
        df_list.append(df_0)
    
    target_average_dict = {}
    
    for df in df_list:
        
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
        
        macs2 = target_average_dict[gene][0]*target_average_dict[gene][1]/len(df_list)
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
    
    source = ['POU5F1']*len(up_target) # TF column
    
    df = pd.DataFrame({'TF': source,
                        'Target Gene': up_target,
                        'macs2': macs2,
                        'q-value': qval})
    df.index += 1
    print(f'Predicted Gene Targets = {len(df.index)}')
    #print(df)
    return(df)

def e7_chip_df_constructor(file):
    """
    1) Takes as argument list objects returned by the chip_reader() functions.
    2) Filters the gene targets.
    3) Returns DataFrames with 4 columns (TF, Target Gene, macs2, q-value).

    Parameters
    ----------
    chip_datatable (list object returned by chip_reader() functions)
    
    Returns
    -------
    pandas DataFrame

    """
    
    import pandas as pd
    
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
    
    df_list = []
    
    for cell_line, experiments in header_dict.items():
        
        header_list = []
        
        for experiment in experiments:
            header = experiment + '|' + cell_line
            header_list.append(header)
        
        header_list.insert(0,'Target_genes')
        
        df = file.loc[:, header_list]
        df['Average'] = df.mean(axis = 1)
        df_0 = df[df['Average']>0]
        df_list.append(df_0)
    
    target_average_dict = {}
    
    for df in df_list:
        
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
        
        macs2 = target_average_dict[gene][0]/target_average_dict[gene][1]
        final_target_average_dict.update({gene:macs2})
   
    final_df = pd.DataFrame(final_target_average_dict.items(), columns = ['Target_genes', 'Average'])
    new_avg = final_df['Average'].tolist()
    tar = final_df['Target_genes'].tolist()
    
    test_list = []
    
    for (x, y) in zip(tar, new_avg):
        test_list.append((x, y))
    
    n = int(1*len(test_list)/100)
    
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
    
    source = ['POU5F1']*len(up_target) # TF column
    
    df = pd.DataFrame({'TF': source,
                      'Target Gene': up_target,
                      'macs2': macs2,
                      'q-value': qval})
    
    df.index += 1
    #print(f'Predicted Gene Targets = {len(df.index)}')
    #print(df)
    return(df)

def deg_df_constructor(datatable):
    """
    1) Takes as an argument list objects returned by deg_reader() function.
    2) Reads gene targets and returns DataFrame with 2 columns (DEG, Expression).

    Parameters
    ----------
    datatable (list object returned by deg_redaer() function)

    Returns
    -------
    pandas DataFrame.

    """
    
    import pandas as pd      
    
    genes = [] # list where DEGs are stored
    high_low = [] # list where expression values are stored (-1,1)
    
    for entry in datatable:
        genes.append(entry[0])
        high_low.append(entry[1])
    
    df = pd.DataFrame({'Genes': genes,
                       'Expression': high_low})
    
    low = len(df[df['Expression']==-1])
    high = len(df[df['Expression']==1])
    
    print(f'Low expressed genes = {low}')
    print(f'High expressed genes = {high}\n')
    return(df)

def chip_df_merger(chipframes):
    """
    1) Takes as argument a list containing all DataFrames returned by
       chip_df_constructor() functions.
    2) Counts the occurance of genes.
    3) Returns merged DataFrames with 4 columns (TF, Target Gene, macs2, q-value).

    Parameters
    ----------
    chipframes : list object which contains all DataFrames returned by
    chip_df_constructor() function.

    Returns
    -------
    df (pandas Dataframe).

    """
    import pandas as pd
    
    all_df = pd.concat(chipframes) # concatenates all dataframes
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
    source = ['POU5F1']*len(gene) # TF column
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
    print(f'Predicted Gene Targets = {len(df.index)}')
    return(df)

def deg_df_merger(degframes):
    """
    1) Takes as argument a list containing all DataFrames returned by
       deg_df_constructor() function.
    2) Counts the occurance of genes.
    3) Returns merged DataFrames with 4 columns (TF, Target Gene, macs2, q-value).

    Parameters
    ----------
    degframes : list object which contains all DataFrames returned by
    deg_df_constructor() function.

    Returns
    -------
    pandas DataFrame

    """
    import pandas as pd
    
    all_df = pd.concat(degframes) # concatenates all dataframes
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
    
    print(f'Integrated DEGs = {len(genes)}')
    print(f'Low expressed genes = {len(lowgenes)}')
    print(f'High expressed genes = {len(highgenes)}\n')
    
    return(df)

def comparer(oct4chip_df, deg_df):
    """
    1) Takes as arguments the data structures returned by chip_df_constructor()
       and deg_df_constructor() function.
    2) Counts the occurance of genes.
    3) Returns merged DataFrames with 4 columns (TF, Target Gene, macs2, q-value).

    Parameters
    ----------
    oct4chip_df (DataFrame returned by chip_df_constructor)
    deg_df (DataFrame returned by deg_df_constructor)

    Returns
    -------
    pandas DataFrame

    """
    import pandas as pd
    
    chip_targets = oct4chip_df['Target Gene'].tolist()
    q_value = oct4chip_df['q-value'].tolist()
    deg_genes = deg_df['DEGs'].tolist()
    expression = deg_df['Expression'].tolist()
    
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
    
    source = ['POU5F1']*len(genes) 
    
    df = pd.DataFrame({'TF': source,
                       'Genes':genes,
                       'q-value':qval,
                       'Expression':exp_1})
    
    print(f'DEGs = {len(commons.keys())} \n')
    return(df)

def e7_comparer(e7chip_df, deg_df):
    """
    1) Takes as arguments the data structures returned by chip_df_constructor()
       and deg_df_constructor() function.
    2) Counts the occurance of genes.
    3) Returns merged DataFrames with 4 columns (TF, Target Gene, macs2, q-value).

    Parameters
    ----------
    e7chip_df (DataFrame returned by chip_df_constructor)
    deg_df (DataFrame returned by deg_df_constructor)

    Returns
    -------
    pandas DataFrame

    """
    import pandas as pd
    
    chip_targets = e7chip_df['Target Gene'].tolist()
    q_value = e7chip_df['q-value'].tolist()
    deg_genes = deg_df['DEGs'].tolist()
    expression = deg_df['Expression'].tolist()
    
    chip_dict = {chip_targets[i]: q_value[i] for i in range(len(chip_targets))}
    deg_dict = {deg_genes[i]: expression[i] for i in range(len(deg_genes))}
    
    commons = {}
    for gene in deg_dict.keys():
        if gene in chip_dict.keys():
            commons.update({gene:[chip_dict[gene],deg_dict[gene]]})
    
    genes = [] # list of TCGA genes that are probably targeted by E7
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
    
    source = ['E7']*len(genes) 
    
    df = pd.DataFrame({'TF': source,
                       'Genes':genes,
                       'q-value':qval,
                       'Expression':exp_1})
    
    print(f'DEGs = {len(commons.keys())} \n')
    return(df)

def gene_comparer(e7_dict, oct4_dict, out_file):
    
    import pandas as pd
    
    common = []
    expression = []
    regulation = []
    
    for gene, exp in e7_dict.items():
        if gene in oct4_dict.keys():
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
    
    print(df)
    #df.to_excel(out_file)
    return(df)
    
def visualizer(dataframe, title):
    """
    1) Takes as argument pandas DataFrame objects.
    2) Establishes a connection with Cytoscape to visualize regulatory networks.

    Parameters
    ----------
    dataframe : pandas DataFrame object returned by comparer() function

    Returns
    -------
    None

    """

    from py2cytoscape.data.cyrest_client import CyRestClient
    from py2cytoscape.cyrest.base import api
    
    # Create an instance of cyREST client.  Default IP is 'localhost', and port number is 1234.
    # cy = CyRestClient() - This default constructor creates connection to http://localhost:1234/v1
    cy = CyRestClient(ip='127.0.0.1', port=1234)
    
    # Cleanup: Delete all existing networks and tables in current Cytoscape session
    # cy.session.delete()
    
    source = dataframe.columns[0]
    target = dataframe.columns[1]
    score = dataframe.columns[2]
    expression = dataframe.columns[3]  
    
    regnet = cy.network.create_from_dataframe(dataframe,
                                              source_col=source,
                                              target_col=target,
                                              interaction_col=expression,
                                              name=title)
    
    # Create a new style
    style1 = cy.style.create('sample_style1')
    
    edge_vps = cy.style.vps.get_edge_visual_props()
    node_vps = cy.style.vps.get_node_visual_props()
    
    arrow_shape = {'High':'Arrow', 'Low':'T'}
    arrow_size = {'High': '10', 'Low':'10'}
    arrow_color = {'High':'green', 'Low':'red'}
    protein_color = {'POU5F1':'purple', 'E7':'purple', 'POU5F1/E7':'purple'}
    basic_settings = {'NODE_FILL_COLOR':'pink', 'NODE_TRANSPARENCY':200}
    style1.update_defaults(basic_settings)
    
    
    style1.create_discrete_mapping(column='interaction',
                                   col_type='String',
                                   vp='EDGE_TARGET_ARROW_SHAPE',
                                   mappings=arrow_shape)
    
    style1.create_discrete_mapping(column='interaction',
                                   col_type='String',
                                   vp='EDGE_TARGET_ARROW_SIZE',
                                   mappings=arrow_size)
    
    style1.create_discrete_mapping(column='interaction',
                                   col_type='String',
                                   vp='EDGE_TARGET_ARROW_UNSELECTED_PAINT',
                                   mappings=arrow_color)
    
    style1.create_discrete_mapping(column='name',
                                   col_type='String',
                                   vp='NODE_FILL_COLOR',
                                   mappings=protein_color)
    
    style1.create_passthrough_mapping(column='name',
                                      col_type='String',
                                      vp='NODE_LABEL')
    
    style1.delete_mapping(vp='EDGE_LABEL')
    
    # Apply the new style
    cy.style.apply(style1, regnet)
    
    # Apply layout
    cy.layout.apply(name='kamada-kawai', network=regnet)
    
    api(namespace="analyzer",command="analyze")

def exporter(dataframe, out_file):
    """
    1) Takes as an argument the merged pandas DataFrame objects.
    2) Outputs a tsv file with 4 columns (Source, Target, Weight, Interaction).

    Parameters
    ----------
    dataframe : pandas DataFrame objects

    Returns
    -------
    tsv file

    """
    import pandas as pd
    
    dataframe.to_excel(out_file)
   

import pandas as pd
import os
import sys

sys.path.append(os.path.realpath('..'))
dirname = os.path.dirname(__file__)

# PPI data files
all_oct4 = os.path.join(dirname, '..\Supplementary File 2\PPI data\All_Oct4.txt')
all_e7 = os.path.join(dirname, '..\Supplementary File 2\PPI data\All_E7.txt')

print('PPI datafiles imported into the algorithm:')
print('1) Integrated Oct4.')
oct4_datatable = ppi_reader(all_oct4)
print('2) Integrated E7.')
e7_datatable = ppi_reader(all_e7)


tss = int(input('Please designate the distance in kb from TSS (1,5 or 10): '))

if tss == 1:
    # Homo sapiens Oct4 ChiP-Atlas gene targets (hg38, 1TSS)
    homo_oct4 = os.path.join(dirname, '..\Supplementary File 1\ChIP-Atlas data\POU5F1\hg38_POU5F1.1.tsv')
    # Mus musculus Oct4 ChiP-Atlas gene targets (mm10, 1TSS)
    mus_oct4 = os.path.join(dirname, '..\Supplementary File 1\ChIP-Atlas data\POU5F1\mm10_Pou5f1.1.tsv')
    
    print('\nChIP-Atlas datafiles imported into the algorithm:')
    print('1) Homo sapiens Oct4 (1kb TSS).')
    hoct4 = oct4_chip_reader(homo_oct4)
    hoct4_df = homo_oct4_chip_df_constructor(hoct4)
    homo_chip_df = hoct4_df
    
    print('2) Mus musculus Oct4 (1kb TSS).')
    moct4 = oct4_chip_reader(mus_oct4)
    moct4_df = mus_oct4_chip_df_constructor(moct4)
    mus_chip_df = moct4_df
    
    chip_df = pd.concat([homo_chip_df, mus_chip_df])
    
    # E7 Homo sapiens top 5 TF ChIP-Atlas gene tagrets (hg38, 1TSS)
    e2f1 = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\hg38_E2F1.1.tsv')
    e2f4 = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\hg38_E2F4.1.tsv')
    erg = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\hg38_ERG.1.tsv')
    myc = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\hg38_MYC.1.tsv')
    tfdp1 = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\hg38_TFDP1.1.tsv')
    
    # E7 Mus musculus top 5 TF ChIP-Atlas gene tagrets (mm10, 1TSS)
    me2f1 = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\mm10_E2f1.1.tsv')
    me2f4 = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\mm10_E2f4.1.tsv')
    merg = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\mm10_Erg.1.tsv')
    mmyc = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\mm10_Myc.1.tsv')
    mtead4 = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\mm10_Tead4.1.tsv')
    
    homo = [e2f1, e2f4, erg, myc, tfdp1]
    mus = [me2f1, me2f4, merg, mmyc, mtead4]
    
    
    e7_chipframes = []
    e7_chipframes_m = []
        
    for file in homo:
        chip = e7_chip_reader(file)
        filtered = e7_chip_df_constructor(chip)
        e7_chipframes.append(filtered)
    
    for file in mus:
        chip = e7_chip_reader(file)
        filtered = e7_chip_df_constructor(chip)
        e7_chipframes_m.append(filtered)
    
    print('3) E7 (Homo sapiens TFs) (1kb TSS).')
    e7_chip_df = chip_df_merger(e7_chipframes)
    print('4) E7 (Mus musculus TFs) (1kb TSS).')
    e7_chip_df_m = chip_df_merger(e7_chipframes_m)
    
    chip_df_b = pd.concat([e7_chip_df, e7_chip_df_m])
    
    # TCGA data files
    c5a2 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-C5-A3HE-01A-21R-A22U-07.txt')
    c5a7 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-C5-A7X8-01A-11R-A36F-07.txt')
    dga2 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-DG-A2KJ-01A-11R-A32Y-07.txt')
    dsa7 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-DS-A7WI-01A-12R-A352-07.txt')
    eka2 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-EK-A2PK-01A-11R-A18M-07.txt')
    eka3gm = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-EK-A3GM-01A-11R-A213-07.txt')
    eka3gn = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-EK-A3GN-01A-11R-A213-07.txt')
    exa1 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-EX-A1H6-01B-11R-A22U-07.txt')
    exa8 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-EX-A8YF-01A-11R-A37O-07.txt')
    hma6 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-HM-A6W2-06A-22R-A33Z-07 .txt')
    ira3 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-IR-A3LI-01A-11R-A32Y-07.txt')
    q1a5 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-Q1-A5R1-01A-11R-A28H-07.txt')
    qia6 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-Q1-A6DV-01A-11R-A32P-07.txt')
    vsa8 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-VS-A8EJ-01A-11R-A37O-07.txt')
    
    tcga_files = [c5a2, c5a7, dga2, dsa7, eka2, eka3gm, eka3gn, exa1, exa8, hma6, ira3, q1a5, qia6, vsa8]
    degframes = []
    
    i = 1
    print('\nTCGA datafiles imported into the algorithm:')
    for file in tcga_files:
        print(f'File No{i}')
        tcga = deg_reader(file)
        frame = deg_df_constructor(tcga)
        degframes.append(frame)
        i += 1
    
    deg_df = deg_df_merger(degframes)
    
    print('Predicted gene targets in the context of Cervical Cancer:')
    print('Homo sapiens Oct4:')
    homo_oct4_cervical = comparer(homo_chip_df, deg_df)
    print('Mus musculus Oct4:')
    mus_oct4_cervical = comparer(mus_chip_df, deg_df)
    oct4_cervical = pd.concat([homo_oct4_cervical, mus_oct4_cervical])
    print('E7 (Homo sapiens TFs)')
    e7_cervical_h = e7_comparer(e7_chip_df, deg_df)
    print('E7 (Mus musculus TFs)')
    e7_cervical_m = e7_comparer(e7_chip_df_m, deg_df)
    e7_cervical = pd.concat([e7_cervical_h, e7_cervical_m])
    
    e7 = e7_cervical['Genes'].tolist()
    
    e7_exp = e7_cervical['Expression'].tolist()
    
    oct4 = oct4_cervical['Genes'].tolist()
    
    oct4_exp = oct4_cervical['Expression'].tolist()
    
    e7_dict = {e7[i]:e7_exp[i] for i in range(len(e7))}
    oct4_dict = {oct4[i]:oct4_exp[i] for i in range(len(oct4))}
    
    print('Common predicted gene targets of Oct4 and E7:')
    oct4_e7_df = gene_comparer(e7_dict, oct4_dict, 'common_1.xlsx')
    l = oct4_e7_df['Target'].tolist()
    
    # print('\nRegulatory Networks visualized in Cytoscape:\n')
    # print('1) Homo sapiens Oct4')
    # visualizer(homo_oct4_cervical, 'Homo sapiens Oct4')
    # print('2) Mus musculus Oct4')
    # visualizer(mus_oct4_cervical, 'Mus musculus Oct4')
    # print('3) Merged Oct4')
    # visualizer(oct4_cervical, 'Merged Oct4')
    # print('4) Homo sapiens - E7')
    # visualizer(e7_cervical_h, 'Homo sapiens - E7')
    # print('5) Mus musculus - E7')
    # visualizer(e7_cervical_m, 'Mus musculus - E7')
    # print('6) Merged E7')
    # visualizer(e7_cervical, 'Merged E7')
    # print('7) Oct4 and E7 Common')
    # visualizer(oct4_e7_df, 'Oct4 and E7 Common')
    
    # exporter(oct4_cervical, 'Supplementary Table 11.xlsx')
    # exporter(e7_cervical, 'Supplementary Table 12.xlsx')
    # exporter(oct4_e7_df, 'Supplementary Table 13.xlsx')

elif tss == 5:
    # Homo sapiens Oct4 ChiP-Atlas gene targets (hg38, 5TSS)
    homo_oct4 = os.path.join(dirname, '..\Supplementary File 1\ChIP-Atlas data\POU5F1\hg38_POU5F1.5.tsv')
    # Mus musculus Oct4 ChiP-Atlas gene targets (mm10, 5TSS)
    mus_oct4 = os.path.join(dirname, '..\Supplementary File 1\ChIP-Atlas data\POU5F1\mm10_Pou5f1.5.tsv')
    
    print('\nChIP-Atlas datafiles imported into the algorithm:')
    print('1) Homo sapiens Oct4 (1kb TSS).')
    hoct4 = oct4_chip_reader(homo_oct4)
    hoct4_df = homo_oct4_chip_df_constructor(hoct4)
    homo_chip_df = hoct4_df
    
    print('2) Mus musculus Oct4 (1kb TSS).')
    moct4 = oct4_chip_reader(mus_oct4)
    moct4_df = mus_oct4_chip_df_constructor(moct4)
    mus_chip_df = moct4_df
    
    chip_df = pd.concat([homo_chip_df, mus_chip_df])
    
    # E7 Homo sapiens top 5 TF ChIP-Atlas gene tagrets (hg38, 5TSS)
    e2f1 = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\hg38_E2F1.5.tsv')
    e2f4 = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\hg38_E2F4.5.tsv')
    erg = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\hg38_ERG.5.tsv')
    myc = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\hg38_MYC.5.tsv')
    tfdp1 = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\hg38_TFDP1.5.tsv')
    
    # E7 Mus musculus top 5 TF ChIP-Atlas gene tagrets (mm10, 5TSS)
    me2f1 = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\mm10_E2f1.5.tsv')
    me2f4 = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\mm10_E2f4.5.tsv')
    merg = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\mm10_Erg.5.tsv')
    mmyc = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\mm10_Myc.5.tsv')
    mtead4 = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\mm10_Tead4.5.tsv')
    
    homo = [e2f1, e2f4, erg, myc, tfdp1]
    mus = [me2f1, me2f4, merg, mmyc, mtead4]
    
    
    e7_chipframes = []
    e7_chipframes_m = []
        
    for file in homo:
        chip = e7_chip_reader(file)
        filtered = e7_chip_df_constructor(chip)
        e7_chipframes.append(filtered)
    
    for file in mus:
        chip = e7_chip_reader(file)
        filtered = e7_chip_df_constructor(chip)
        e7_chipframes_m.append(filtered)
    
    print('3) E7 (Homo sapiens TFs) (1kb TSS).')
    e7_chip_df = chip_df_merger(e7_chipframes)
    print('4) E7 (Mus musculus TFs) (1kb TSS).')
    e7_chip_df_m = chip_df_merger(e7_chipframes_m)
    
    chip_df_b = pd.concat([e7_chip_df, e7_chip_df_m])
    
    # TCGA data files
    c5a2 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-C5-A3HE-01A-21R-A22U-07.txt')
    c5a7 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-C5-A7X8-01A-11R-A36F-07.txt')
    dga2 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-DG-A2KJ-01A-11R-A32Y-07.txt')
    dsa7 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-DS-A7WI-01A-12R-A352-07.txt')
    eka2 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-EK-A2PK-01A-11R-A18M-07.txt')
    eka3gm = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-EK-A3GM-01A-11R-A213-07.txt')
    eka3gn = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-EK-A3GN-01A-11R-A213-07.txt')
    exa1 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-EX-A1H6-01B-11R-A22U-07.txt')
    exa8 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-EX-A8YF-01A-11R-A37O-07.txt')
    hma6 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-HM-A6W2-06A-22R-A33Z-07 .txt')
    ira3 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-IR-A3LI-01A-11R-A32Y-07.txt')
    q1a5 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-Q1-A5R1-01A-11R-A28H-07.txt')
    qia6 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-Q1-A6DV-01A-11R-A32P-07.txt')
    vsa8 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-VS-A8EJ-01A-11R-A37O-07.txt')
    
    tcga_files = [c5a2, c5a7, dga2, dsa7, eka2, eka3gm, eka3gn, exa1, exa8, hma6, ira3, q1a5, qia6, vsa8]
    degframes = []
    
    i = 1
    print('\nTCGA datafiles imported into the algorithm:')
    for file in tcga_files:
        print(f'File No{i}')
        tcga = deg_reader(file)
        frame = deg_df_constructor(tcga)
        degframes.append(frame)
        i += 1
    
    deg_df = deg_df_merger(degframes)
    
    print('Predicted gene targets in the context of Cervical Cancer:')
    print('Homo sapiens Oct4:')
    homo_oct4_cervical = comparer(homo_chip_df, deg_df)
    print('Mus musculus Oct4:')
    mus_oct4_cervical = comparer(mus_chip_df, deg_df)
    oct4_cervical = pd.concat([homo_oct4_cervical, mus_oct4_cervical])
    print('E7 (Homo sapiens TFs)')
    e7_cervical_h = e7_comparer(e7_chip_df, deg_df)
    print('E7 (Mus musculus TFs)')
    e7_cervical_m = e7_comparer(e7_chip_df_m, deg_df)
    e7_cervical = pd.concat([e7_cervical_h, e7_cervical_m])
    
    e7 = e7_cervical['Genes'].tolist()
    
    e7_exp = e7_cervical['Expression'].tolist()
    
    oct4 = oct4_cervical['Genes'].tolist()
    
    oct4_exp = oct4_cervical['Expression'].tolist()
    
    e7_dict = {e7[i]:e7_exp[i] for i in range(len(e7))}
    oct4_dict = {oct4[i]:oct4_exp[i] for i in range(len(oct4))}
    
    print('Common predicted gene targets of Oct4 and E7:')
    oct4_e7_df = gene_comparer(e7_dict, oct4_dict, 'common_5.xlsx')
    l = oct4_e7_df['Target'].tolist()
    
    # print('\nRegulatory Networks visualized in Cytoscape:\n')
    # print('1) Homo sapiens Oct4')
    # visualizer(homo_oct4_cervical, 'Homo sapiens Oct4')
    # print('2) Mus musculus Oct4')
    # visualizer(mus_oct4_cervical, 'Mus musculus Oct4')
    # print('3) Merged Oct4')
    # visualizer(oct4_cervical, 'Merged Oct4')
    # print('4) Homo sapiens - E7')
    # visualizer(e7_cervical_h, 'Homo sapiens - E7')
    # print('5) Mus musculus - E7')
    # visualizer(e7_cervical_m, 'Mus musculus - E7')
    # print('6) Merged E7')
    # visualizer(e7_cervical, 'Merged E7')
    # print('7) Oct4 and E7 Common')
    # visualizer(oct4_e7_df, 'Oct4 and E7 Common')
    
    # exporter(oct4_cervical, 'Supplementary Table 11.xlsx')
    # exporter(e7_cervical, 'Supplementary Table 12.xlsx')
    # exporter(oct4_e7_df, 'Supplementary Table 13.xlsx')

elif tss == 10:
    # Homo sapiens Oct4 ChiP-Atlas gene targets (hg38, 10TSS)
    homo_oct4 = os.path.join(dirname, '..\Supplementary File 1\ChIP-Atlas data\POU5F1\hg38_POU5F1.10.tsv')
    # Mus musculus Oct4 ChiP-Atlas gene targets (mm10, 10TSS)
    mus_oct4 = os.path.join(dirname, '..\Supplementary File 1\ChIP-Atlas data\POU5F1\mm10_Pou5f1.10.tsv')
    
    print('\nChIP-Atlas datafiles imported into the algorithm:')
    print('1) Homo sapiens Oct4 (10kb TSS).')
    hoct4 = oct4_chip_reader(homo_oct4)
    hoct4_df = homo_oct4_chip_df_constructor(hoct4)
    homo_chip_df = hoct4_df
    
    print('2) Mus musculus Oct4 (10kb TSS).')
    moct4 = oct4_chip_reader(mus_oct4)
    moct4_df = mus_oct4_chip_df_constructor(moct4)
    mus_chip_df = moct4_df
    
    chip_df = pd.concat([homo_chip_df, mus_chip_df])
    chip_df.reset_index(inplace = True)
    # print(chip_df)
    
    # E7 Homo sapiens top 5 TF ChIP-Atlas gene tagrets (hg38, 10TSS)
    e2f1 = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\hg38_E2F1.10.tsv')
    e2f4 = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\hg38_E2F4.10.tsv')
    erg = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\hg38_ERG.10.tsv')
    myc = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\hg38_MYC.10.tsv')
    tfdp1 = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\hg38_TFDP1.10.tsv')
    
    # E7 Mus musculus top 5 TF ChIP-Atlas gene tagrets (mm10, 10TSS)
    me2f1 = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\mm10_E2f1.10.tsv')
    me2f4 = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\mm10_E2f4.10.tsv')
    merg = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\mm10_Erg.10.tsv')
    mmyc = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\mm10_Myc.10.tsv')
    mtead4 = os.path.join(dirname,'..\Supplementary File 1\ChIP-Atlas data\E7\mm10_Tead4.10.tsv')
    
    homo = [e2f1, e2f4, erg, myc, tfdp1]
    mus = [me2f1, me2f4, merg, mmyc, mtead4]
    
    
    e7_chipframes = []
    e7_chipframes_m = []
        
    for file in homo:
        chip = e7_chip_reader(file)
        filtered = e7_chip_df_constructor(chip)
        e7_chipframes.append(filtered)
    
    for file in mus:
        chip = e7_chip_reader(file)
        filtered = e7_chip_df_constructor(chip)
        e7_chipframes_m.append(filtered)
    
    print('3) E7 (Homo sapiens TFs) (10kb TSS).')
    e7_chip_df = chip_df_merger(e7_chipframes)
    print('4) E7 (Mus musculus TFs) (10kb TSS).')
    e7_chip_df_m = chip_df_merger(e7_chipframes_m)
    
    chip_df_b = pd.concat([e7_chip_df, e7_chip_df_m])
    chip_df_b.reset_index(inplace = True)
    # print(chip_df_b)
    
    # TCGA data files
    c5a2 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-C5-A3HE-01A-21R-A22U-07.txt')
    c5a7 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-C5-A7X8-01A-11R-A36F-07.txt')
    dga2 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-DG-A2KJ-01A-11R-A32Y-07.txt')
    dsa7 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-DS-A7WI-01A-12R-A352-07.txt')
    eka2 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-EK-A2PK-01A-11R-A18M-07.txt')
    eka3gm = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-EK-A3GM-01A-11R-A213-07.txt')
    eka3gn = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-EK-A3GN-01A-11R-A213-07.txt')
    exa1 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-EX-A1H6-01B-11R-A22U-07.txt')
    exa8 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-EX-A8YF-01A-11R-A37O-07.txt')
    hma6 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-HM-A6W2-06A-22R-A33Z-07 .txt')
    ira3 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-IR-A3LI-01A-11R-A32Y-07.txt')
    q1a5 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-Q1-A5R1-01A-11R-A28H-07.txt')
    qia6 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-Q1-A6DV-01A-11R-A32P-07.txt')
    vsa8 = os.path.join(dirname,'..\Supplementary File 1\TCGA data\CESC_TCGA-VS-A8EJ-01A-11R-A37O-07.txt')
    
    tcga_files = [c5a2, c5a7, dga2, dsa7, eka2, eka3gm, eka3gn, exa1, exa8, hma6, ira3, q1a5, qia6, vsa8]
    degframes = []
    
    i = 1
    print('\nTCGA datafiles imported into the algorithm:')
    for file in tcga_files:
        print(f'File No{i}')
        tcga = deg_reader(file)
        frame = deg_df_constructor(tcga)
        degframes.append(frame)
        i += 1
    
    deg_df = deg_df_merger(degframes)
    
    print('Predicted gene targets in the context of Cervical Cancer:')
    print('Homo sapiens Oct4:')
    homo_oct4_cervical = comparer(homo_chip_df, deg_df)
    print('Mus musculus Oct4:')
    mus_oct4_cervical = comparer(mus_chip_df, deg_df)
    oct4_cervical = pd.concat([homo_oct4_cervical, mus_oct4_cervical])
    oct4_cervical.reset_index(0)
    print('E7 (Homo sapiens TFs)')
    e7_cervical_h = e7_comparer(e7_chip_df, deg_df)
    print('E7 (Mus musculus TFs)')
    e7_cervical_m = e7_comparer(e7_chip_df_m, deg_df)
    e7_cervical = pd.concat([e7_cervical_h, e7_cervical_m])
    e7_cervical.reset_index(0)
    
    e7 = e7_cervical['Genes'].tolist()
    
    e7_exp = e7_cervical['Expression'].tolist()
    
    oct4 = oct4_cervical['Genes'].tolist()
    
    oct4_exp = oct4_cervical['Expression'].tolist()
    
    e7_dict = {e7[i]:e7_exp[i] for i in range(len(e7))}
    oct4_dict = {oct4[i]:oct4_exp[i] for i in range(len(oct4))}
    
    print('Common predicted gene targets of Oct4 and E7:')
    oct4_e7_df = gene_comparer(e7_dict, oct4_dict, 'common_10.xlsx')
    l = oct4_e7_df['Target'].tolist()
    
    print('\nRegulatory Networks visualized in Cytoscape:\n')
    # print('1) Homo sapiens Oct4')
    # visualizer(homo_oct4_cervical, 'Homo sapiens Oct4')
    # print('2) Mus musculus Oct4')
    # visualizer(mus_oct4_cervical, 'Mus musculus Oct4')
    print('3) Merged Oct4')
    visualizer(oct4_cervical, 'Merged Oct4')
    # print('4) Homo sapiens - E7')
    # visualizer(e7_cervical_h, 'Homo sapiens - E7')
    # print('5) Mus musculus - E7')
    # visualizer(e7_cervical_m, 'Mus musculus - E7')
    print('6) Merged E7')
    visualizer(e7_cervical, 'Merged E7')
    # print('7) Oct4 and E7 Common')
    # visualizer(oct4_e7_df, 'Oct4 and E7 Common')
    
    exporter(oct4_cervical, 'Supplementary Table 11.xlsx')
    exporter(e7_cervical, 'Supplementary Table 12.xlsx')
    exporter(oct4_e7_df, 'Supplementary Table 13.xlsx')
