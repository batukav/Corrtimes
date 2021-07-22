# Corrtimes
Project to investigate the dynamics of current lipid force fields

This project is initiated by H. Antila and M. Miettinen but we welcome anyone interested to participate. For publications arising from this project, the tentative author order is H. Antila first and M. Miettinen last, while others are placed in between. In case of significant contributions this is reconsidered in a fair manner.

Manuscript *Using open data to rapidly benchmark biomolecular simulations: Phospholipid conformational dynamics* has been accepted for publication in Journal of Chemical Information and Modelling, [preprint at bioRxiv](https://www.biorxiv.org/content/10.1101/2020.11.09.374850v2).

The databank port contains the files necessary for running the correlation times within the NMRLipids databank format. The biggest change to the original script is the removal of Gromacs rotacf function. Instead, now we are using an equivalent Numpy implementation.
