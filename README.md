# IC50 prediction

Predicting IC50 and Uncovering More Effective Drugs through Machine Learning Models

## Overview

This jupyter notebook creates and uses a machine learning model to predict the IC50 and search for drugs with a lower IC50.

![Flowchart!](https://github.com/ivanbarraborrella/ic50prediction/blob/main/diagram_2.png)

First, the target will be selected. For this the ChEMBL Database will be used.

Next we use the lipinski descriptor to explore the data and see if there is any significant difference between the active and inactive molecules. 

In this step we will select a number of molecules, in this case 10, and we will create analogues to these, specifically in the code nitrogen atoms are changed for carbon atoms, but this can be changed and create analogues in many other ways.

Now we will use the PaDEL descriptors software. PaDEL-Descriptor is like a translator that converts the structure of a molecule into a language that computers can understand and process.

Imagine you have a molecule, which is a complex three-dimensional structure. For a human, it can be difficult to understand all the properties of this molecule just by looking at it. But if we use PaDEL-Descriptor, we can get a list of numbers (called descriptors) that represent different aspects of the molecule, such as its size, shape, and the presence of certain groups of atoms.

These descriptors can then be used to predict properties of the molecule, such as its reactivity, its toxicity, or whether it could be a good drug for a certain disease. This is very useful in fields like medicinal chemistry and bioinformatics, where efficient ways to analyze and compare molecules are needed.

Now, we create the X and Y matrix to create the machine learning model. We compare the different algorithms and choose the one that best suits us. In the next step we create and save our machine learning model.



## References

[1] Chanin Nantasenamat. https://www.youtube.com/dataprofessor

[2] Yap, C.W. (2011). PaDEL-Descriptor: An open source software to calculate molecular descriptors and fingerprints. Journal of Computational Chemistry, 32(7), 1466-1474. doi: 10.1002/jcc.21707
