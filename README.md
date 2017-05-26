# graphmolsim
Set of scripts enabling to use maximum common substructure and graph edit distance in ligand-based virtual screening.

This set of scripts enables to test different graph algorithms in ligand-based virtual screening.
The main approaches are Maximum Common Subgraph (MCS) and Graph Edit Distance (GED).
Both approaches can be run through the script main.sh.

Requirements
--------
In order to make the scripts working, the RDKit chemoinformatics library and path to it has to be exported
in the RDBASE environment variable.
The scripts are written for Python v. 2.7. Python v. 3 has not yet been tested.
Java 1.7 has to be present on the system. 
This software has been primarily written for Linux, the starting scripts are written in Bash.
However, the bash scripts can be replaced and rest of the software is platform independent.

Datasets
--------
The input data can be downloaded from these addresses:
http://siret.ms.mff.cuni.cz/skoda/datasets/8.0-8.5.zip

http://siret.ms.mff.cuni.cz/skoda/datasets/8.5-9.0.zip

http://siret.ms.mff.cuni.cz/skoda/datasets/9.0-9.5.zip

http://siret.ms.mff.cuni.cz/skoda/datasets/9.8-1.0.zip
The .zip files have to be saved and upacked in the same directory as the scripts.

MCS
--------
The fmcs function counting MCS directly from RDKit is used.

GED
--------
Graph Matching Toolkit (GMT) from Kaspar Risen has been used for counting GED.
More information on http://www.fhnw.ch/wirtschaft/iwi/gmt
