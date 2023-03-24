### Bachelor Thesis
### Computational Analysis of viral-host protein interactions: focus on E7 & POU5F1.

Researched and implemented a Python pipeline for the analysis of in-house and publicly available genomics and proteomics data. Regulatory and protein interaction networks were constructed and visualized using Plotly and NetworkX packages to extract information about a protein complexes.

*CytoNet.ipynb & CytoNet.py*
- Read the input PPI datafiles and returns an appropriate data structure.
- Calculate an interaction score based on IntAct Database scoring method.
- Extract the common interactors of the input datafiles.
- Visualize the PPI networks using plotly and networkx.

*ComplexFinder.ipynb & ComplexFinder.py*
- Read the PPI datafiles and return appropriate data structures.
- Read the HGNC datafiles and return appropriate datastructures.
- Compare data and visualize the results as bar charts using matplotlib.

*ChIPNet.ipynb & ChIPNet.py*
- Read the PPI, ChIP-Atlas and TCGA input datafiles.
- Return an appropriate data structure for each input file.
- Compare data and return the predicted gene targets.

