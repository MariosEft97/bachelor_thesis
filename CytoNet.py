import pandas as pd

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
                if score < 0.166:
                    interaction_indication.append('Weak')
                if score >= 0.166 and score < 0.332:
                    interaction_indication.append('Medium')
                if score >= 0.498:
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
