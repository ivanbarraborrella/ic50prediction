# IC50 prediction

Predicting IC50 and Uncovering More Effective Drugs through Machine Learning Models

## Overview

This jupyter notebook creates and uses a machine learning model to predict the IC50 and search for drugs with a lower IC50.

![Flowchart!](https://github.com/ivanbarraborrella/ic50prediction/blob/main/diagram_2.png)

Initially, we select our target using the ChEMBL Database and prepare the necessary data. 

Subsequently, we employ the Lipinski descriptor to explore the data, aiming to identify any significant differences between the active and inactive molecules. 

In the following step, we select a set of molecules, in this case, ten, and create analogues for these. Specifically, the code changes nitrogen atoms to carbon atoms in these molecules. However, this can be modified to create analogues in various other ways.

We then utilize the PaDEL-Descriptor software. PaDEL-Descriptor serves as a translator, converting the structure of a molecule into a language that computers can comprehend and process.

Consider a molecule as a complex three-dimensional structure. For a human, understanding all the properties of this molecule merely by looking at it can be challenging. However, with PaDEL-Descriptor, we can obtain a list of numbers, known as descriptors, representing different aspects of the molecule, such as its size, shape, and the presence of certain groups of atoms.

These descriptors can subsequently be used to predict various properties of the molecule, such as its reactivity, toxicity, or potential as a drug for a specific disease. This is particularly useful in fields like medicinal chemistry and bioinformatics, where efficient methods to analyze and compare molecules are required.

Next, we construct the X and Y matrices to create a machine learning model. We compare various algorithms and select the one that best fits our needs. Following this, we create and save our machine learning model.

We then predict the pIC50 for our analogues and select the promising ones. 

Finally, we perform several feedback cycles, creating analogues of our analogues, until we are left with the final candidates.

## References

[1] Chanin Nantasenamat. https://www.youtube.com/dataprofessor  https://github.com/dataprofessor

[2] Yap, C.W. (2011). PaDEL-Descriptor: An open source software to calculate molecular descriptors and fingerprints. Journal of Computational Chemistry, 32(7), 1466-1474. doi: 10.1002/jcc.21707
