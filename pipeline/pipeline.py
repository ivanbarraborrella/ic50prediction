import pandas as pd
from chembl_webresource_client.new_client import new_client
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

import data_collection
import exploratory_data

def pipeline(protein):
    
    data_collection.data_collection(protein)

    exploratory_data.exploratory_data1('bioactivity_data_curated.csv')


pipeline('sos')
