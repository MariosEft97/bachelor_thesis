# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 13:35:50 2021

@author: Marios
"""

"""
CytoNet_2.3.0

The current version of the algorithm does the following:
    
    1) Reads the input PPI datafiles and returns an appropriate data structure.
    2) Calculates an interaction score based on IntAct Database scoring method.
    3) Extracts the common interactors of the input datafiles.
    4) Automatically visualizes the PPI networks in Cytoscape platform.
    5) Utilizes stringApp to perform functional enrichment analysis.
    6) Utilizes stringApp to expand the networks.
    7) Outputs results in tsv files.

    *** Cytoscape platform (3.8.2) needs to be open for the algorithm to run.
    *** In each algorithm run only one visualizer function can be called.
    *** Request limitations permit only one automatic enrichment visualization
        per algorithm run.
    
    Algorithm Runs:
    
    1) Calls visualizer() function to visualize:
        a) Homo sapiens Oct4 interactome.
        b) HPV16 E7 interactome.
        c) Integrated Oct4 interactome.')
        d) Integrated E7 interactome.
    
    2) Calls visualizer_1() function to visualize:
        a) Homo sapiens Oct4 & HPV16 E7 interactome.
        b) Integrated Oct4 & E7 interactome.
    
    3) Calls visualizer_2() function to visualize:
        a) Homo sapiens Oct4 & HPV16 E7 Common Interactors.
        b) Integrated Oct4 & E7 Common Interactors.
      
       Calls fea() to perform Functional Enrichment Analysis
   
    4) Calls expansion() function to expand:
        a) Integrated Oct4 & E7 Common Interactors.    
    
"""

def reader(datafile):
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
    
    print(f'PPIs = {interactions}.\n')
    return(datatable)
    
def calculator(datatable):
    """
    1) Takes as argument the list object returned by reader() function.
    2) Calculates an interaction score for each entry of the datatable based on
       IntAct scoring methodology.
    3) Creates non-redundant pandas DataFrames.

    Parameters
    ----------
    datatable : list object returned by reader() function.

    Returns
    -------
    dataframe (pandas DataFrame object)

    """
    import pandas as pd
    
    # Data structures used by the algorithm for the interaction score calculation.
    biochemical = ['pull down (MI:0096)',
                   'anti tag coimmunoprecipitation (MI:0007)',
                   'ubiquitinase assay (MI:0997)',
                   'peptide array (MI:0081)',
                   'proximity ligation assay (MI:0813)',
                   'protein array (MI:0089)',
                   'filamentous phage display (MI:0048)'
                   'biochemical (MI:0401)',
                   'enzymatic study (MI:0415)',
                   'competition binding (MI:0405)',
                   'coimmunoprecipitation (MI:0019)',
                   'cosedimentation (MI:0027)',
                   'tandem affinity purification (MI:0676)',
                   'cross-linking study (MI:0030)', 
                   'anti bait coimmunoprecipitation (MI:0006)',
                   'affinity chromatography technology (MI:0004)', 
                   'chromatin immunoprecipitation assay (MI:0402)'
                   'comigration in non denaturing gel electrophoresis (MI:0404)']
    
    biophysical = ['molecular sieving (MI:0071)',
                   'fluorescent resonance energy transfer (MI:0055)',
                   'x-ray crystallography (MI:0114)',
                   'fluorescence technology (MI:0051)',
                   'detection by mass spectrometry (MI:0943)']
    
    pca = ['two hybrid pooling approach (MI:0398)',
           'two hybrid (MI:0018)',
           'two hybrid array (MI:0397)', 
           'beta lactamase complementation (MI:0011)']
    
    imaging_tech = ['imaging technique (MI:0428)',
                    'fluorescence microscopy (MI:0416)',
                    'confocal microscopy (MI:0663)']
    
    computational = ['predictive text mining (MI:0087)']
    
    methodtypescores = {'association (MI:0914)':1,
                        'physical association (MI:0915)':2,
                        'direct interaction (MI:0407)':5,
                        'ubiquitination reaction (MI:0220)':5,
                        'phosphorylation reaction (MI:0220)':5, 
                        'colocalization (MI:0403)':0.2,
                        'predicted interaction (MI:1110)':0.5}
    
    scorelist = [] # list object where initial ppi scores are saved
    interactions = {} # key = pair of interactors, value = summed ppi score
    prot_occurance = {} # key = interactor, value = occurance
    publications = {} # list with publications for each interaction
    methods = {} # list with methods for each interaction
    prot_publication = {} # key = interactor, value = No of publications
    prot_methods = {} # key = interactor, value = No of methods
    id_list = {} # key = interactor protein & value = UniProt ID
    tppi = 0 # total protein - protein interactions in the input file
    
    for entry in datatable:
        tppi += 1
        
        entrynameA = entry[0].split('_') # Protein entry name of interactor A
        speciesA = entrynameA[1] # Species of interactor A
        uniprotA_ID = entry[1]  # UniProt ID of interactor A
        geneA = entry[2] # HGNC gene name of interactor A
        
        entrynameB = entry[3].split('_') # Protein entry name of interactor B
        speciesB = entrynameB[1] # Species of interactor B
        uniprotB_ID = entry[4] # UniProt ID of interactor B
        geneB = entry[5] # HGNC gene name of interactor B
        
        method_type = entry[6] # MICV method type of interaction
        method = entry[7] # MICV method of interaction
        publication = entry[8] # Publication (Authors & PMID)
        source = entry[9] # Primary source of interaction
        
        pair = (geneA, geneB) # Pairs of gene names of interacting proteins
        
        # UniProt IDs of interactor proteins are saved in a dictionary.
        if geneB in id_list.keys():
            pass
        else:
            id_list.update({geneB:uniprotB_ID})
        
        if geneB not in prot_occurance.keys():                
            # key = interactor, value = number of occurance of interactor
            prot_occurance.update({geneB:1})            
        else:
            prot_occurance[geneB] += 1
        
        if geneB not in methods.keys():
            methods.update({geneB: [method]})
        else:
            if method not in methods[geneB]:
                methods[geneB].append(method)
               
        if geneB not in publications.keys():
            publications.update({geneB: [publication]})
        else:
            if publication not in publications[geneB]:
                publications[geneB].append(publication)
        
        # key = interaction pair, value = number of distinct methods   
        prot_methods.update({geneB:len(methods[geneB])})
        # key = interaction pair, value = number of distinct publications
        prot_publication.update({geneB:len(publications[geneB])})
        
        # Score calculation of Oct4 interactions.
        if geneA == 'POU5F1':            
            # variables for interaction-score calculation
            mtscore = 0 # method type score
            mscore = 0  # method score
            score = 0  # interaction score
            
            # Interaction score calculation procedure:
            # 1) IntAct method type and method scores are summed to calculate the interaction score.
            # 2) Final score is calculated by summing the scores of identical interactions.
            if entry[6] in methodtypescores.keys():
                mtscore = methodtypescores[entry[6]]
                
            if entry[7] in biochemical:
                mscore += 3
            elif entry[7] in biophysical:
                mscore += 3
            elif entry[7] in pca:
                mscore += 2
            elif entry[7] in imaging_tech:
                mscore += 0.6
            elif entry[7] in computational:
                mscore += 0.5
            else:
                mscore += 0.1
            
            # For Oct4 interactions observed in other species rathen than Homo sapiens
            # score is multiplied by the identity percentage of proteins as determined by BlastP
            # Mus musculus Oct4 identity = 0.83
            # Xenopus laevis Oct4 identity =  0.56
            if speciesA == 'MOUSE':
                score = (mtscore + mscore) * 0.83
            elif speciesA == 'XENLA':
                score = (mtscore + mscore) * 0.56
            else:
                score = mtscore + mscore
            
            # Each interaction score is stored in a list.
            scorelist.append(score)
            
            # Final score of non-redundant interactions is calculated.
            if pair in interactions:
                interactions[pair] += score
            else:
                interactions.update({pair:score})
        
        elif geneA == 'E7':                
            # variables for interaction-score calculation
            mtscore = 0 # method type score
            mscore = 0  # method score
            score = 0  # interaction score
            
            # Interaction score calculation procedure:
            # 1) IntAct method type and method scores are summed to calculate the interaction score.
            # 2) Final score is calculated by summing the scores of identical interactions.
            if entry[6] in methodtypescores.keys():
                mtscore = methodtypescores[entry[6]]
                
            if entry[7] in biochemical:
                mscore += 3
            elif entry[7] in biophysical:
                mscore += 3
            elif entry[7] in pca:
                mscore += 2
            elif entry[7] in imaging_tech:
                mscore += 0.6
            elif entry[7] in computational:
                mscore += 0.5
            else:
                mscore += 0.1
            
            # For E7 interactions observed in other strains rather than HPV16
            # score is multiplied by the identity percentage of proteins as determined by BlastP
            if speciesA == 'HPV1':
                score = (mtscore + mscore) * 0.38
            elif speciesA == 'HPV2A':
                score = (mtscore + mscore) * 0.46
            elif speciesA == 'HPV04':
                score = (mtscore + mscore) * 0.30
            elif speciesA == 'HPV05':
                score = (mtscore + mscore) * 0.35
            elif speciesA == 'HPV6B':
                score = (mtscore + mscore) * 0.57
            elif speciesA == 'HPV08':
                score = (mtscore + mscore) * 0.37
            elif speciesA == 'HPV11':
                score = (mtscore + mscore) * 0.56
            elif speciesA == 'HPV17':
                score = (mtscore + mscore) * 0.41
            elif speciesA == 'HPV18':
                score = (mtscore + mscore) * 0.42
            elif speciesA == 'HPV20':
                score = (mtscore + mscore) * 0.32
            elif speciesA == 'HPV25':
                score = (mtscore + mscore) * 0.36
            elif speciesA == 'HPV31':
                score = (mtscore + mscore) * 0.74
            elif speciesA == 'HPV33':
                score = (mtscore + mscore) * 0.61
            elif speciesA == 'HPV37':
                score = (mtscore + mscore) * 0.38
            elif speciesA == 'HPV38':
                score = (mtscore + mscore) * 0.33
            elif speciesA == 'HPV45':
                score = (mtscore + mscore) * 0.47
            elif speciesA == 'HPV55':
                score = (mtscore + mscore) * 0.47
            elif speciesA == 'HPV57':
                score = (mtscore + mscore) * 0.46
            elif speciesA == 'HPV74':
                score = (mtscore + mscore) * 0.52
            elif speciesA == 'HPV76':
                score = (mtscore + mscore) * 0.1
            elif speciesA == 'HPV92':
                score = (mtscore + mscore) * 0.1
            elif speciesA == 'HPV98':
                score = (mtscore + mscore) * 0.1
            elif speciesA == 'HPV':
                score = (mtscore + mscore) * 0.1
            else:
                score = mtscore + mscore
            
            # Each interaction score is stored in a list.
            scorelist.append(score)
            
            # Final score of non redundant interactions is calculated.
            if pair in interactions:
                interactions[pair] += score
            else:
                interactions.update({pair:score})
                
    pairs = interactions.keys()
    pairscores = interactions.values()
    
    # final scores are rounded to 2 decimal places
    roundedscores = []
    for score in pairscores:
        roundedscores.append(round(score, 2))
    
    # final scores are summed with the numbers of distinct publications 
    finalpairscores = []
    for pair,score in interactions.items():
        if pair[1] in prot_occurance.keys():
            score = score + prot_publication[pair[1]] # + prot_methods[pair[1]]
            finalpairscores.append(score)
    
    # final scores are normalized to 0-1 range (Min-Max methodology)
    normalized = []
    maximum = max(finalpairscores)
    minimum = min(finalpairscores)
    for score in finalpairscores:
        score = (score-minimum)/(maximum-minimum)
        normalized.append(round(float(score),3))
    
    interactorAlist = [] # from each pair interactor A is saved in this list object.
    interactorBlist = [] # from each pair interactor B is saved in this list object.
    interaction_indication = [] # For each interaction an indication (Weak, Medium, Strong) is saved.
    id_listA = [] # UniProt IDs of interactors A.
    id_listB = [] # UniProt IDs of interactors B.
    methodsAB = [] # Number of distinct methods.
    publicationsAB = []  # Number of distinct publications.
    
    # Interactors are separated and saved in different list objects.
    for pair in pairs:
        interactorA = pair[0]
        interactorAlist.append(interactorA)
        interactorB = pair[1]
        interactorBlist.append(interactorB)
    
    for gene in interactorBlist:
        if gene in id_list.keys():
            id_listB.append(id_list[gene])
        if gene in prot_methods.keys():
            methodsAB.append(prot_methods[gene])
        if gene in prot_publication.keys():
            publicationsAB.append(prot_publication[gene])
    
    id_listA = ['Q01860']*len(id_listB)
    
    # For each interaction score an affinity indication is devised.
    for score in normalized:
        if score < 0.333:
            interaction_indication.append('Weak')
        if score >= 0.333 and score < 0.666:
            interaction_indication.append('Medium')
        if score >= 0.666:
            interaction_indication.append('Strong')
    
    dataframe = pd.DataFrame({'Source': interactorAlist,
                              'Source ID': id_listA,
                              'Target': interactorBlist,
                              'Target ID': id_listB,
                              'Weight': normalized,
                              'Interaction': interaction_indication,
                              'Methods': methodsAB,
                              'Publications': publicationsAB})
    dataframe.index += 1
    print(f'Protein interactors = {len(dataframe.index)}')
    
    return(dataframe)

def common(oct4_dataframe, e7_dataframe):
    """
    1) Takes as arguments the DataFrame objects returned by calculator().
    2) Compares them and returns new DataFrame with their common entities.

    Parameters
    ----------
    oct4_dataframe : pandas DataFrame object returned by calculator()
    e7_dataframe : pandas DataFrame object returned by calculator()

    Returns
    -------
    dataframe (pandas DataFrame object)

    """
    import pandas as pd
    
    common = []
    
    oct4list = oct4_dataframe['Target'].tolist() # oct4 ppi targets saved in a list
    e7list = e7_dataframe['Target'].tolist() # e7 ppi targets saved in a list
    
    # Oct4 and E7 ppi targets are compared.
    for prot in oct4list:
        if prot in e7list:
            common.append(prot)
    
    # Pandas Dataframe for visualization purposes is constructed.
    df = oct4_dataframe[oct4_dataframe.Target.isin(e7_dataframe.Target)]
    df1 = e7_dataframe[e7_dataframe.Target.isin(oct4_dataframe.Target)]
    dataframe = pd.concat([df,df1])
    
    dataframe.index += 1
    print(f'Protein interactors = {len(dataframe.index)}')
    
    return(dataframe)

def visualizer(dataframe, title):
    """
    1) Takes as argument DataFrame objects returned from calculator.
    2) Establishes a connection with Cytoscape platform.
    3) Visualizes the PPI networks.
    
    *** visualizer() function is used for Oct4 and E7 interactomes.

    Parameters
    ----------
    dataframe : pandas DataFrame object returned by calculator()

    Returns
    -------
    None

    """

    from py2cytoscape.data.cyrest_client import CyRestClient
    from py2cytoscape.cyrest.base import api
    
    # Create an instance of cyREST client.
    # Default IP is 'localhost', and port number is 1234.
    # The default constructor creates connection to http://localhost:1234/v1
    cy = CyRestClient(ip='127.0.0.1', port=1234)
    
    # Cleanup: Delete all existing networks and tables in current Cytoscape session
    # cy.session.delete()
    
    source = dataframe.columns[0]
    target = dataframe.columns[2]
    score = dataframe.columns[4]
    interaction = dataframe.columns[5]
        
    ppin = cy.network.create_from_dataframe(dataframe,
                                            source_col=source,
                                            target_col=target,
                                            interaction_col=interaction,
                                            name=title)
    
    # Grab list of nodes and edges
    node_suids = ppin.get_nodes()
    edge_suids = ppin.get_edges()
    
    
    # Create a new style
    style1 = cy.style.create('sample_style1')
    
    edge_vps = cy.style.vps.get_edge_visual_props()
    node_vps = cy.style.vps.get_node_visual_props()
    # print(edge_vps)
    # print(node_vps)
    
    # Discrete mapping: Simply prepare key-value pairs and send it
    edge_color = {'Weak': 'purple', 'Medium': 'navy', 'Strong': 'black'}
    edge_width = {'Weak': '0.5', 'Medium': '0.75', 'Strong': '1'}
    
    basic_settings = {'NODE_FILL_COLOR':'red', 'NODE_TRANSPARENCY':200}
    style1.update_defaults(basic_settings)
    
    # Apply the new style
    cy.style.apply(style1, ppin)
    
    # Discrete mappings
    style1.create_discrete_mapping(column='interaction',
                                   col_type='String',
                                   vp='EDGE_STROKE_UNSELECTED_PAINT',
                                   mappings=edge_color)
    
    style1.create_discrete_mapping(column='interaction', 
                                   col_type='String',
                                   vp='EDGE_WIDTH',
                                   mappings=edge_width)
    
    # Passthrough mapping
    style1.create_passthrough_mapping(column='name',
                                      col_type='String',
                                      vp='NODE_LABEL')
    
    # Apply layout
    cy.layout.apply(name='kamada-kawai', network=ppin)
    
    api(namespace="analyzer",command="analyze")

def visualizer_1(dataframe, title):
    """
    1) Takes as argument DataFrame objects returned from calculator.
    2) Establishes a connection with Cytoscape platform.
    3) Visualizes the PPI networks.
    
    *** visualizer_1() function is used for merged Oct4 and E7 interactomes.

    Parameters
    ----------
    dataframe : pandas DataFrame object returned by calculator()

    Returns
    -------
    None

    """

    from py2cytoscape.data.cyrest_client import CyRestClient
    from py2cytoscape.cyrest.base import api
    
    # Create an instance of cyREST client.
    # Default IP is 'localhost', and port number is 1234.
    # The default constructor creates connection to http://localhost:1234/v1
    cy = CyRestClient(ip='127.0.0.1', port=1234)
    
    # Cleanup: Delete all existing networks and tables in current Cytoscape session
    # cy.session.delete()
    
    source = dataframe.columns[0]
    target = dataframe.columns[2]
    score = dataframe.columns[4]
    interaction = dataframe.columns[5]
    
    ppin = cy.network.create_from_dataframe(dataframe,
                                            source_col=source,
                                            target_col=target,
                                            interaction_col=interaction,
                                            name=title)
    
    # Grab list of nodes and edges
    node_suids = ppin.get_nodes()
    edge_suids = ppin.get_edges()
    
    # Create a new style
    style1 = cy.style.create('sample_style1')
    
    edge_vps = cy.style.vps.get_edge_visual_props()
    node_vps = cy.style.vps.get_node_visual_props()
    # print(edge_vps)
    # print(node_vps)
    
    # Discrete mapping: Simply prepare key-value pairs and send it
    edge_color = {'Weak': 'purple', 'Medium': 'navy', 'Strong': 'black'}
    edge_width = {'Weak': '1', 'Medium': '1', 'Strong': '1'}
    protein_color_1 = {'POU5F1':'yellow','E7':'yellow'}
    protein_size = {'POU5F1':'200','E7':'100' }
    node_label_size = {'POU5F1':'50','E7':'40' }
    
    basic_settings = {'NODE_FILL_COLOR':'purple'}
    style1.update_defaults(basic_settings)
    
    # Apply the new style
    cy.style.apply(style1, ppin)
    
    # Discrete mappings
    style1.create_discrete_mapping(column='interaction',
                                   col_type='String',
                                   vp='EDGE_STROKE_UNSELECTED_PAINT',
                                   mappings=edge_color)
    
    style1.create_discrete_mapping(column='interaction',
                                   col_type='String',
                                   vp='EDGE_WIDTH',
                                   mappings=edge_width)
    
    style1.create_discrete_mapping(column='name', 
                                   col_type='String',
                                   vp='NODE_FILL_COLOR',
                                   mappings=protein_color_1)
    
    style1.create_discrete_mapping(column='name',
                                   col_type='String',
                                   vp='NODE_SIZE',
                                   mappings=protein_size)
    
    style1.create_discrete_mapping(column='name',
                                   col_type='String',
                                   vp='NODE_LABEL_FONT_SIZE',
                                   mappings=node_label_size)
    
    # Passthrough mapping
    style1.create_passthrough_mapping(column='name',
                                      col_type='String',
                                      vp='NODE_LABEL')
    
    # Apply layout
    cy.layout.apply(name='kamada-kawai', network=ppin)
    
    api(namespace="analyzer",command="analyze")
    
def visualizer_2(dataframe, title):
    """
    1) Takes as argument DataFrame objects returned from calculator.
    2) Establishes a connection with Cytoscape platform.
    3) Visualizes the PPI networks.
    
    *** visualizer_2() function is used for Oct4 and E7 common interactors.

    Parameters
    ----------
    dataframe : pandas DataFrame object returned by calculator()

    Returns
    -------
    None

    """

    from py2cytoscape.data.cyrest_client import CyRestClient
    from py2cytoscape.cyrest.base import api
    
    # Create an instance of cyREST client.
    # Default IP is 'localhost', and port number is 1234.
    # The default constructor creates connection to http://localhost:1234/v1
    cy = CyRestClient()
    
    # Cleanup: Delete all existing networks and tables in current Cytoscape session
    # cy.session.delete()
    
    source = dataframe.columns[0]
    target = dataframe.columns[2]
    score = dataframe.columns[4]
    expression = dataframe.columns[5]
        
    ppin = cy.network.create_from_dataframe(dataframe,
                                            source_col=source,
                                            target_col=target,
                                            interaction_col=expression,
                                            name=title)
    
    # Grab list of nodes and edges
    node_suids = ppin.get_nodes()
    edge_suids = ppin.get_edges()
    
    # Create a new style
    style1 = cy.style.create('sample_style1')
    
    edge_vps = cy.style.vps.get_edge_visual_props()
    node_vps = cy.style.vps.get_node_visual_props()
    # print(edge_vps)
    # print(node_vps)
    
    # Discrete mapping: Simply prepare key-value pairs and send it
    edge_color = {'Weak': 'purple', 'Medium': 'navy', 'Strong': 'black'}
    edge_width = {'Weak': '1', 'Medium': '1', 'Strong': '1'}
    protein_color = {'POU5F1':'yellow','E7':'yellow' }
    
    basic_settings = {'NODE_FILL_COLOR':'red', 'NODE_TRANSPARENCY':200}
    style1.update_defaults(basic_settings)
    
    # Apply the new style
    cy.style.apply(style1, ppin)
    
    # Discrete mapping
    style1.create_discrete_mapping(column='interaction',
                                   col_type='String',
                                   vp='EDGE_STROKE_UNSELECTED_PAINT',
                                   mappings=edge_color)
    
    style1.create_discrete_mapping(column='interaction',
                                   col_type='String',
                                   vp='EDGE_WIDTH',
                                   mappings=edge_width)
    
    style1.create_discrete_mapping(column='name',
                                   col_type='String',
                                   vp='NODE_FILL_COLOR',
                                   mappings=protein_color)
    
    # Passthrough mapping
    style1.create_passthrough_mapping(column='name',
                                      col_type='String',
                                      vp='NODE_LABEL')
    
    # Apply layout
    cy.layout.apply(name='cose', network=ppin)
    
    api(namespace="analyzer",command="analyze")
    
def fea(dataframe, term):
    """
    1) Takes as argument the DataFrames returned by calculator() function.
    2) Utilizes stringApp to perform Functional Enrichment Analysis.

    Parameters
    ----------
    dataframe : pandas DataFrame object returned by calculator()

    Returns
    -------
    None.

    """
    from py2cytoscape.cyrest.base import api
    import time
    
    ci = dataframe['Target'].tolist()
    ci_unique = []
    
    for prot in ci:
        if prot == 'E7':
            pass
        else:
            if prot in ci_unique:
                pass
            else:
                ci_unique.append(prot)
    
    ci_query = ','.join(ci_unique)
    #print(ci_query)
    
    
    api(namespace="string",
        command="protein query",
        PARAMS={"query":ci_query,
                "cutoff":"0.1",
                "limit":"0",
                "includesViruses": "true",
                "network": "current",
                "newNetName": "Common Interactors Network",
                "species":"Homo sapiens",
                "taxonID": "9606"})
    
    api(namespace='string',
        command='add nodes',
        PARAMS={"query": "E7",
                "cutoff": "0.1",
                "includesViruses": "true",
                "limit": "1",
                "network": "current",
                "species": "Human papillomavirus type 16",
                "taxonID": "337041"})
    
    api(namespace='network',
        command='delete',
        PARAMS={'network': 'CURRENT',
                'nodeList': '337041.VE6_HPV16'})
    
    api(namespace='layout',
        command='kamada-kawai')

    
    api(namespace='string',
        command='retrieve enrichment',
        PARAMS={"allNetSpecies": "Homo sapiens",
                "background": "genome",
                "selectedNodesOnly": "false",
                "view": "string"})
    
    api(namespace='string',
        command='show enrichment')
    
    api(namespace='string',
        command='filter enrichment',
        PARAMS={"categories": term,
                "overlapCutoff": "0.5",
                "removeOverlapping": "true"})
    
    api(namespace='string',command='show charts')
    
    api(namespace="analyzer",command="analyze")

def expansion(dataframe, num):
    """
    1) Takes as argument the DataFrames returned by calculator() function.
    2) Utilizes stringApp to expand the given PPI network.

    Parameters
    ----------
    dataframe : pandas DataFrame object returned by calculator()

    Returns
    -------
    None.

    """
    from py2cytoscape.cyrest.base import api
    import time
    
    ci = dataframe['Target'].tolist()
    ci_unique = []
    
    for prot in ci:
        if prot == 'E7':
            pass
        else:
            if prot in ci_unique:
                pass
            else:
                ci_unique.append(prot)
    
    ci_query = ','.join(ci_unique)
    #print(ci_query)
    
    
    api(namespace="string",
        command="protein query",
        PARAMS={"query":ci_query,
                "cutoff":"0.1",
                "limit":"0",
                "includesViruses": "true",
                "network": "current",
                "newNetName": "Common Interactors Network",
                "species":"Homo sapiens",
                "taxonID": "9606"})

    api(namespace='string',
        command='add nodes',
        PARAMS={"query": "E7",
                "cutoff": "0.1",
                "includesViruses": "true",
                "limit": "1",
                "network": "current",
                "species": "Human papillomavirus type 16",
                "taxonID": "337041"})
    
    api(namespace='network',
        command='delete',
        PARAMS={'network': 'CURRENT',
                'nodeList': '337041.VE6_HPV16'})
    
    api(namespace='string',
        command='expand',
        PARAMS={"additionalNodes": num,
                "network": "current",
                "nodeTypes": 'Homo sapiens',
                "selectivityAlpha": "0.7"})
    
    
    api(namespace='layout',command='cose')
    
    api(namespace="analyzer",command="analyze")
    
def exporter(dataframe, out_file):
    """
    1) Takes as argument DataFrame objects.
    2) Outputs a tsv file with 8 columns.

    Parameters
    ----------
    dataframe : pandas DataFrame objects returned by calculator()
    out_file: name of the output file

    Returns
    -------
    None

    """
   #  import csv
        
   # # A tsv data file with interactors is created.
   #  with open(out_file,'w', newline = '') as output_file:
   #      tsvfile = csv.writer(output_file, delimiter = '\t')
   #      tsvfile.writerow(['Source', 'Source_ID', 'Target', 'Target_ID', 'Weight', 'Interaction', 'Methods', 'Publications'])  
        
   #      for entry in dataframe.iloc:
   #          tsvfile.writerow([entry[0], entry[1], entry[2], entry[3], entry[4], entry[5], entry[6], entry[7]]) 
   
    import pandas as pd
    dataframe.to_excel(out_file)
   
   
import pandas as pd 
import os
import sys

sys.path.append(os.path.realpath('..'))
dirname = os.path.dirname(__file__)

all_oct4 = os.path.join(dirname, '..\Supplementary File 2\PPI data\All_Oct4.txt')
all_e7 = os.path.join(dirname, '..\Supplementary File 2\PPI data\All_E7.txt')
hs_oct4 = os.path.join(dirname, '..\Supplementary File 2\PPI data\Oct4_Homo_sapiens.txt')
hpv16_e7 = os.path.join(dirname, '..\Supplementary File 2\PPI data\E7_HPV16.txt')


print('PPI datafiles imported into the algorithm:')
print('1) Homo sapiens Oct4.')
hsoct4_datatable = reader(hs_oct4)
print('2) HPV16 E7.')
hpv16e7_datatable = reader(hpv16_e7)
print('3) Integrated Oct4.')
oct4_datatable = reader(all_oct4)
print('4) Integrated E7.')
e7_datatable = reader(all_e7)

print('\nNon-redundant PPI dataframes constructed:')

print('\n1) Homo sapiens Oct4.')
hsoct4_dataframe = calculator(hsoct4_datatable)

print('\n2) HPV16 E7.')
hpv16e7_dataframe = calculator(hpv16e7_datatable)

print('\n3) Integrated Oct4.')
oct4_dataframe = calculator(oct4_datatable)

print('\n4) Integrated E7.')
e7_dataframe = calculator(e7_datatable)

print('\n5) Homo sapiens Oct4 and HPV16 E7 Common Interactors.')
hs_hpv16_ci = common(hsoct4_dataframe, hpv16e7_dataframe)
print('\n6) Integrated Oct4 and E7 Common Interactors.')
ci = common(oct4_dataframe, e7_dataframe)

dataframe = pd.concat([hsoct4_dataframe,hpv16e7_dataframe])
dataframe1 = pd.concat([oct4_dataframe, e7_dataframe])

# In each algorithm run, only one visualizer function can be called.
algorithm_run = int(input('Please designate the algorithm run (1,2,3 or 4): '))

if algorithm_run == 1:
    # visualizer function call
    
    from py2cytoscape.data.cyrest_client import CyRestClient
    
    # Create an instance of cyREST client.
    # Default IP is 'localhost', and port number is 1234.
    # The default constructor creates connection to http://localhost:1234/v1
    cy = CyRestClient(ip='127.0.0.1', port=1234)
    
    # Cleanup: Delete all existing networks and tables in current Cytoscape session
    cy.session.delete()
    
    print('\nPPI networks visualized in Cytoscape platform:')
    print('1) Homo sapiens Oct4 interactome.')
    visualizer(hsoct4_dataframe, 'Homo sapiens Oct4 interactome')
    print('2) HPV16 E7 interactome.')
    visualizer(hpv16e7_dataframe, 'HPV16 E7 interactome')
    print('3) Integrated Oct4 interactome.')
    visualizer(oct4_dataframe, ' Integrated Oct4 interactome')
    print('4) Integrated E7 interactome.')
    visualizer(e7_dataframe, 'Integrated E7 interactome')

elif algorithm_run == 2:
    # visualizer_1 function call
    
    from py2cytoscape.data.cyrest_client import CyRestClient
    
    # Create an instance of cyREST client.
    # Default IP is 'localhost', and port number is 1234.
    # The default constructor creates connection to http://localhost:1234/v1
    cy = CyRestClient(ip='127.0.0.1', port=1234)
    
    # Cleanup: Delete all existing networks and tables in current Cytoscape session
    cy.session.delete()
    
    print('1) Homo sapiens Oct4 & HPV16 E7 interactome.')
    visualizer_1(dataframe, 'Homo sapiens Oct4 & HPV16 E7 interactome')
    print('2) Integrated Oct4 & E7 interactome.')
    visualizer_1(dataframe1, 'Integrated Oct4 & E7 interactome')

elif algorithm_run == 3:
    # visualizer_2 function call + Functional Enrichment Analysis
        
    from py2cytoscape.data.cyrest_client import CyRestClient
    
    # Create an instance of cyREST client.
    # Default IP is 'localhost', and port number is 1234.
    # The default constructor creates connection to http://localhost:1234/v1
    cy = CyRestClient(ip='127.0.0.1', port=1234)
    
    # Cleanup: Delete all existing networks and tables in current Cytoscape session
    cy.session.delete()
    
    print('1) Homo sapiens Oct4 & HPV16 E7 Common Interactors.')
    visualizer_2(hs_hpv16_ci, 'Homo sapiens Oct4 & HPV16 E7 Common Interactors')
    print('2) Integrated Oct4 & E7 Common Interactors.')
    visualizer_2(ci, 'Integrated Oct4 & E7 Common Interactors')

    # Request limitations permit only one automatic enrichment visualization per run.
    # The desired terms can be manually selected in each algorithm run.
    print('\nFunctional Enrichment Analysis of Oct4 & E7 Common Interactors Network:')
    print('\nGO Process = 1, GO Function = 2, GO Component = 3, KEGG Pathways = 4')
    analysis = input('Please type 1,2,3 or 4: ')
    print('\nPerforming Analysis...')
    fea(ci, analysis)
    print('\nAnalysis performed!')
    
elif algorithm_run == 4:
    # Expansion call
    
    from py2cytoscape.data.cyrest_client import CyRestClient
    
    # Create an instance of cyREST client.
    # Default IP is 'localhost', and port number is 1234.
    # The default constructor creates connection to http://localhost:1234/v1
    cy = CyRestClient(ip='127.0.0.1', port=1234)
    
    # Cleanup: Delete all existing networks and tables in current Cytoscape session
    cy.session.delete()
    
    # PPI network expansions
    exp = ['20', '30', '40', '50']
    
    
    print('\nExpansion of Oct4 & E7 Common Interactors Network:')
    for num in exp:
        print(f'Network is expanded by {num} new interactors.')
        expansion(ci, num)

# print('\nExcel output files of the algorithm:')
# print('1) Supplementary Table 1')
# exporter(oct4_dataframe, 'Supplementary Table 1.xlsx')
# print('2) Supplementary Table 2')
# exporter(e7_dataframe, 'Supplementary Table 2.xlsx')
# print('3) Supplementary Table 3')
# exporter(ci, 'Supplementary Table 3.xlsx')
