### Bachelor Thesis
### Computational Analysis of viral-host protein interactions: focus on E7 & POU5F1.

Researched and implemented a Python pipeline for the analysis of in-house and publicly available genomics and proteomics data. Regulatory and protein interaction networks were constructed and visualized within CytoScape to extract information about a specific protein complex.

CytoNet.py functionality:
- Reads the input PPI datafiles and returns an appropriate data structure.
- Calculates an interaction score based on IntAct Database scoring method.
- Extracts the common interactors of the input datafiles.
- Automatically visualizes the PPI networks in Cytoscape platform.
- Utilizes stringApp to perform functional enrichment analysis.
- Utilizes stringApp to expand the networks.
- Outputs results in tsv files.

ComplexFinder.py functionality:
- Reads the PPI datafiles and returns an appropriate data structure.
- Reads the HGNC datafiles and returns an appropriate datastructure.
- Compares them and visualizes the results as bar charts.

ChIPNet.py functionality:
- Reads the PPI, ChIP-Atlas and TCGA input datafiles.
- Returns an appropriate data structure for each input file.
- Compares them and returns the predicted gene targets.
- Visualizes results in Cytoscape platform.

