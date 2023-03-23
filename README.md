### Bachelor Thesis
### Computational Analysis of viral-host protein interactions: focus on E7 & POU5F1.

Researched and implemented a Python pipeline for the analysis of in-house and publicly available genomics and proteomics data. Regulatory and protein interaction networks were constructed and visualized using Plotly and NetworkX packages to extract information about a specific protein complex.

*CytoNet.py*
- Reads the input PPI datafiles and returns an appropriate data structure.
- Calculates an interaction score based on IntAct Database scoring method.
- Extracts the common interactors of the input datafiles.
- Visualizes the PPI networks using plotly and networkx.

*ComplexFinder.py*
- Reads the PPI datafiles and returns an appropriate data structure.
- Reads the HGNC datafiles and returns an appropriate datastructure.
- Compares them and visualizes the results as bar charts.

*ChIPNet.py*
- Reads the PPI, ChIP-Atlas and TCGA input datafiles.
- Returns an appropriate data structure for each input file.
- Compares them and returns the predicted gene targets.

