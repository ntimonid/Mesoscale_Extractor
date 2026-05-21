# Mesoscale_Extractor
Given a source and a target mouse brain area for which the user wishes to delineate projection patterns, Mesoscale Extractor downloads publicly available neuronal morphologies, whose soma is located in the source and whose axonal terminals are located in the target areas, and it estimates meso-scale connectivity statistics similar to how it is computed in tract-tracing experiments such as in Oh et al 2014.  

## References
This code has been used in the following published work:  
Timonidis, Nestor, et al. "Translating single-neuron axonal reconstructions into meso-scale connectivity statistics in the mouse somatosensory thalamus." Frontiers in neuroinformatics 17 (2023): 1272243. doi: https://doi.org/10.3389/fninf.2023.1272243

## Description 
Files:  
mesoscale_extractor_example.ipynb: examplar Jupyter Notebook that illustrates the utility of the tool
mesoscales_extractor.py: the main code behind the tool
cfg.py: a configuration script containing all the libraries used by the tool

Directories:   

libraries: a number of auxiliary Python scripts that are being used as back-end libraries for the tool.
  1. convertAllenSpace.py: contains a number of functions related to the affine transformation operations needed to register neuronal morphologies to a reference space
  2. NeuronMorphology.py: class used for storing, manipulating and analysing neuronal morphologies as Python class objects
  3. utils.py: a collection of back-end functions that are important for the tool.  
     
     
atlas_files: a number of nrrd, pkl and json files containing pre-computed information that is fundamental for the spatial orientation of the data to the Common Coordinate Framework developed by the Allen Institute, collection of morphologies that have been grouped based on the common anatomical location of their soma area which dramatically speeds up the search for neurons needed to extract mesoscale statistics and anatomical information related to the distinct individual barrels of the mouse primary somatosensory whiskers cortex (wS1). 


