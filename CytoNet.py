import pandas as pd
import plotly.graph_objects as go
import networkx as nx
import matplotlib.pyplot as plt
from itertools import zip_longest
import plotly.io as pio
pio.renderers.default = 'colab'

# ----------------------------------------------------------------------------------------------------

def reader(data_filepath: str) -> list:
    
    """
    The function reads the tsv input data files and stores each line as a list in a list object.
    
    Parameters:
        data_filepath (str): filepath to data
    
    Returns: 
        datatable (list): data structure to save data
    """
    
    if not isinstance(data_filepath, str):
        raise TypeError("data_filepath must be specified as a string")
    
    else:
    
        interactions = 0 # variable to count interactions in each data file

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

def calculator(datatable: list, protein: str) -> pd.DataFrame:
    
    """
    The function takes as input the list object returned by reader() function. It calculates an
    interaction score for each entry of the datatable based on IntAct scoring methodology and returns
    a non-redundant pandas DataFrame with the results.

    Parameters:
        datatable (list): data structure returned from reader() function
        protein (str): string denonting protein of interest (options: oct4 or e7)

    Returns:
        dataframe (pd.DataFrame): data structure containing interaction score results

    """
    
    if not isinstance(datatable, list):
        raise TypeError("datatable must be specified as a list")
    
    elif not isinstance(protein, str):
        raise TypeError("protein must be specified as a string\noptions: oct4 or e7")
    
    else:
    
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
                # print(prot_publication[pair[1]])
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

        if protein == 'oct4':
            id_listA = ['Q01860']*len(id_listB)

            # For each interaction score an affinity indication is devised.
            for score in normalized:
                if score < 0.333:
                    interaction_indication.append('Weak')
                if score >= 0.333 and score < 0.666:
                    interaction_indication.append('Medium')
                if score >= 0.666:
                    interaction_indication.append('Strong')
        else:
            id_listA = ['P03129']*len(id_listB)

            # For each interaction score an affinity indication is devised.
            for score in normalized:
                if score < 0.333:
                    interaction_indication.append('Weak')
                if score >= 0.333 and score < 0.666:
                    interaction_indication.append('Medium')
                if score >= 0.666:
                    interaction_indication.append('Strong')


        df = pd.DataFrame(
            {
                'Source': interactorAlist,
                'Source ID': id_listA,
                'Target': interactorBlist,
                'Target ID': id_listB,
                'Weight': normalized,
                'Interaction': interaction_indication,
                'Methods': methodsAB,
                'Publications': publicationsAB
            }
        )
        
        df.index += 1

        return(df)

# ----------------------------------------------------------------------------------------------------

def common_interactors(protein_a_df: pd.DataFrame, protein_b_df: pd.DataFrame) -> pd.DataFrame:
    
    """
    The function takes as arguments the DataFrame objects returned from calculator() function.
    It compares the two input DataFrames and returns new DataFrame containing only the common interactors.

    Parameters:
        protein_a_df (pd.DataFrame): data structure returned from calculator() function
        protein_b_df (pd.DataFrame): data structure returned from calculator() function

    Returns:
        df (pd.DataFrame): data structure containing only common interactors

    """
    
    if not isinstance(protein_a_df, pd.DataFrame):
        raise TypeError("protein_a_df must be specified as a pandas DataFrame")
    
    elif not isinstance(protein_b_df, pd.DataFrame):
        raise TypeError("protein_b_df must be specified as a pandas DataFrame")
    
    else:
    
        common_interactors = []

        oct4_interactors_list = protein_a_df['Target'].tolist() # oct4 ppi targets saved in a list
        e7_interactors_list = protein_b_df['Target'].tolist() # e7 ppi targets saved in a list

        # Oct4 and E7 ppi targets are compared.
        for interactor in oct4_interactors_list:
            if interactor in e7_interactors_list:
                common_interactors.append(interactor)

        # Pandas Dataframe for visualization purposes is constructed.
        df_a = protein_a_df[protein_a_df.Target.isin(protein_b_df.Target)]
        df_b = protein_b_df[protein_b_df.Target.isin(protein_a_df.Target)]
        df = pd.concat([df_a, df_b])

        df.index += 1

        return(df)
    
# ----------------------------------------------------------------------------------------------------

def common_interactor_network(df: pd.DataFrame) -> None:
    
    """
    The function takes as arguments the DataFrame objects returned from common_interactors() function
    and visualizes the network using plotly and networkx.

    Parameters:
        df (pd.DataFrame): data structure returned from common_interactors() function

    Returns:
        None

    """
    
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be specified as a pandas DataFrame")
    
    else:
        
        source_nodes = list(df["Source"])
        target_nodes = list(df["Target"])
        unique_nodes = pd.Series(source_nodes+target_nodes).drop_duplicates().tolist()
        weights = list(df["Weight"])
        interactions = list(df["Interaction"])
    
        G = nx.Graph()
        
        # add network nodes
        G.add_nodes_from(unique_nodes)
        
        # add network weighted edges
        weighted_edges = [(source, target, weight) for source, target, weight in zip(source_nodes, target_nodes, weights)]
        G.add_weighted_edges_from(weighted_edges)
        
        # positions for all nodes - seed for reproducibility
        pos = nx.spring_layout(G, seed=0)
        
        edge_x = []
        edge_y = []

        for edge in G.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.append(x0)
            edge_x.append(x1)
            edge_x.append(None)
            edge_y.append(y0)
            edge_y.append(y1)
            edge_y.append(None)

        edge_trace = go.Scatter(
            x=edge_x, y=edge_y,
            line=dict(width=0.5),
            hoverinfo='text',
            mode='lines'
        )

        node_x = []
        node_y = []
        for node in G.nodes():
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)

        node_trace = go.Scatter(
            x=node_x,
            y=node_y,
            mode='markers',
            hoverinfo='text',
            marker=dict(showscale=True, colorscale='YlGnBu', reversescale=True, color=[], size=20,
                        colorbar=dict(thickness=15, title='Protein Interactors', xanchor='left', titleside='right'), line_width=2))
        
        node_adjacencies = []
        node_text = []


        for node, adjacencies in zip(G.nodes, G.adjacency()):
            node_adjacencies.append(len(adjacencies[1]))
            # node_text.append(node+", "+str(len(adjacencies[1]))+" interactors")

        pou5f1_interactions = {}
        e7_interactions = {}

        for edge in G.edges(data=True):
            if edge[0] == "POU5F1" and edge[1] == "E7":
                e7_interactions.update({edge[0]: edge[2]["weight"]})

        for edge in G.edges(data=True):
            if edge[0] == "POU5F1":
                pou5f1_interactions.update({edge[1]: edge[2]["weight"]})
            elif edge[0] == "E7":
                e7_interactions.update({edge[1]: edge[2]["weight"]})

        sorted_pou5f1_interactions = sorted(pou5f1_interactions.items(), key=lambda pair: list(G.nodes).index(pair[0]))
        sorted_e7_interactions = sorted(e7_interactions.items(), key=lambda pair: list(G.nodes).index(pair[0]))

        for node, pou5f1_node, e7_node in zip(list(G.nodes()), sorted_pou5f1_interactions, sorted_e7_interactions):
            node_text.append(f"{pou5f1_node[0]} --- POU5F1: {str(pou5f1_node[1])}<br>{e7_node[0]} --- E7: {str(e7_node[1])}")

        node_trace.marker.color = node_adjacencies
        node_trace.text = node_text
        
        edge_dash = ""
        edge_text = []
        for edge in (G.edges(data=True)):
            edge_text.append(edge[2]["weight"])
            if 0 <= edge[2]["weight"] < 0.33:
                edge_dash += "5px "
            elif 0.33 <= edge[2]["weight"] < 0.66:
                edge_dash += "5px "
            elif edge[2]["weight"] >= 0.66:
                edge_dash += "5px "

        # edge_trace.line.dash = edge_dash.rstrip()
        edge_trace.text = edge_text
        
        
        fig = go.Figure(data=[edge_trace, node_trace],
             layout=go.Layout(
                title='<br>POU5F1 & E7 common protein interactors',
                titlefont_size=16,
                showlegend=False,
                hovermode='closest',
                margin=dict(b=20,l=5,r=5,t=40),
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                )
        
        fig.show(renderer="colab")
        
        return(None)
    
# ----------------------------------------------------------------------------------------------------

def interaction_network(protein_a_df: pd.DataFrame, protein_b_df: pd.DataFrame, common_interactors_df: pd.DataFrame) -> None:
    
    """
    The function takes as arguments the DataFrame objects returned from calculator() function and visualizes
    the interaction network using plotly and networkx.

    Parameters:
        protein_a_df (pd.DataFrame): data structure returned from calculator() function
        protein_b_df (pd.DataFrame): data structure returned from calculator() function
        common_interactors_df (pd.DataFrame): data structure returned from common_interactors() function

    Returns:
        None

    """
    
    if not isinstance(protein_a_df, pd.DataFrame):
        raise TypeError("protein_a_df must be specified as a pandas DataFrame")
    
    elif not isinstance(protein_b_df, pd.DataFrame):
        raise TypeError("protein_b_df must be specified as a pandas DataFrame")
    
    elif not isinstance(common_interactors_df, pd.DataFrame):
        raise TypeError("common_interactors_df must be specified as a pandas DataFrame")
    
    else:
        
        concatenated_df = pd.concat([protein_a_df, protein_b_df], axis=0)
        
        source_nodes = list(concatenated_df["Source"])
        target_nodes = list(concatenated_df["Target"])
        unique_nodes = pd.Series(source_nodes+target_nodes).drop_duplicates().tolist()
        weights = list(concatenated_df["Weight"])
        interactions = list(concatenated_df["Interaction"])
        oct4_targets = list(protein_a_df["Target"].unique())
        e7_targets = list(protein_b_df["Target"].unique())
        common_interactors = pd.Series(common_interactors_df["Target"]).drop_duplicates().tolist()

        G = nx.Graph()

        # add network nodes
        G.add_nodes_from(unique_nodes)

        # add network weighted edges
        weighted_edges = [(source, target, weight) for source, target, weight in zip(source_nodes, target_nodes, weights)]
        G.add_weighted_edges_from(weighted_edges)

        # positions for all nodes - seed for reproducibility
        pos = nx.random_layout(G, seed=0)

        edge_x = []
        edge_y = []

        for edge in G.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.append(x0)
            edge_x.append(x1)
            edge_x.append(None)
            edge_y.append(y0)
            edge_y.append(y1)
            edge_y.append(None)

        edge_trace = go.Scatter(
            x=edge_x,
            y=edge_y,
            name="Interactions",
            line=dict(width=0.5, color="#888"),
            hoverinfo='text',
            mode='lines'
        )

        node_x = []
        node_y = []
        for node in G.nodes():
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)

        node_trace = go.Scatter(
            x=node_x,
            y=node_y,
            name="POU5F1 / E7",
            mode='markers',
            hoverinfo='text',
            marker=dict(color=[], size=20))

        node_adjacencies = []
        node_text = []
        node_color = []


        for node, adjacencies in zip(G.nodes, G.adjacency()):
            node_adjacencies.append(len(adjacencies[1]))
            node_text.append(node+", "+str(len(adjacencies[1]))+" interactors")

            if node == "POU5F1" or node == "E7":
                node_color.append("#ffffd9")
            elif node in common_interactors:
                node_color.append("#7fcdbb")
            elif node in oct4_targets:
                node_color.append("#1d91c0")
            elif node in e7_targets:
                node_color.append("#081d58")

        # aliceblue, antiquewhite, aqua, aquamarine, azure,
        # beige, bisque, black, blanchedalmond, blue,
        # blueviolet, brown, burlywood, cadetblue,
        # chartreuse, chocolate, coral, cornflowerblue,
        # cornsilk, crimson, cyan, darkblue, darkcyan,
        # darkgoldenrod, darkgray, darkgrey, darkgreen,
        # darkkhaki, darkmagenta, darkolivegreen, darkorange,
        # darkorchid, darkred, darksalmon, darkseagreen,
        # darkslateblue, darkslategray, darkslategrey,
        # darkturquoise, darkviolet, deeppink, deepskyblue,
        # dimgray, dimgrey, dodgerblue, firebrick,
        # floralwhite, forestgreen, fuchsia, gainsboro,
        # ghostwhite, gold, goldenrod, gray, grey, green,
        # greenyellow, honeydew, hotpink, indianred, indigo,
        # ivory, khaki, lavender, lavenderblush, lawngreen,
        # lemonchiffon, lightblue, lightcoral, lightcyan,
        # lightgoldenrodyellow, lightgray, lightgrey,
        # lightgreen, lightpink, lightsalmon, lightseagreen,
        # lightskyblue, lightslategray, lightslategrey,
        # lightsteelblue, lightyellow, lime, limegreen,
        # linen, magenta, maroon, mediumaquamarine,
        # mediumblue, mediumorchid, mediumpurple,
        # mediumseagreen, mediumslateblue, mediumspringgreen,
        # mediumturquoise, mediumvioletred, midnightblue,
        # mintcream, mistyrose, moccasin, navajowhite, navy,
        # oldlace, olive, olivedrab, orange, orangered,
        # orchid, palegoldenrod, palegreen, paleturquoise,
        # palevioletred, papayawhip, peachpuff, peru, pink,
        # plum, powderblue, purple, red, rosybrown,
        # royalblue, rebeccapurple, saddlebrown, salmon,
        # sandybrown, seagreen, seashell, sienna, silver,
        # skyblue, slateblue, slategray, slategrey, snow,
        # springgreen, steelblue, tan, teal, thistle, tomato,
        # turquoise, violet, wheat, white, whitesmoke,
        # yellow, yellowgreen

        node_trace.marker.color = node_color
        node_trace.text = node_text

        edge_dash = ""
        edge_text = []
        for edge in (G.edges(data=True)):
            edge_text.append(edge[2]["weight"])
            if 0 <= edge[2]["weight"] < 0.33:
                edge_dash += "5px "
            elif 0.33 <= edge[2]["weight"] < 0.66:
                edge_dash += "5px "
            elif edge[2]["weight"] >= 0.66:
                edge_dash += "5px "

        # edge_trace.line.dash = edge_dash.rstrip()
        edge_trace.text = edge_text

        ci = go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            name="Common Interactors",
            marker=dict(size=20, color="#7fcdbb", symbol='circle'),
        )

        pou5f1 = go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            name="POU5F1 Interactors",
            marker=dict(size=20, color="#1d91c0", symbol='circle'),
        )

        e7 = go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            name="E7 Interactors",
            marker=dict(size=20, color="#081d58", symbol='circle'),
        )


        fig = go.Figure(data=[edge_trace, node_trace, ci, pou5f1, e7],
             layout=go.Layout(
                title='<br>POU5F1 & E7 PPI network',
                titlefont_size=16,
                showlegend=True,
                hovermode='closest',
                margin=dict(b=20,l=5,r=5,t=40),
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                )

        fig.show(renderer="colab")
        
        return(None)
    
# ----------------------------------------------------------------------------------------------------