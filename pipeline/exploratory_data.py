import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

# Inspired by: https://codeocean.com/explore/capsules?query=tag:data-curation

def lipinski(smiles, verbose=False):

    moldata= []
    for elem in smiles:
        mol=Chem.MolFromSmiles(elem) 
        moldata.append(mol)
       
    baseData= np.arange(1,1)
    i=0  
    for mol in moldata:        
       
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)
           
        row = np.array([desc_MolWt,
                        desc_MolLogP,
                        desc_NumHDonors,
                        desc_NumHAcceptors])   
    
        if(i==0):
            baseData=row
        else:
            baseData=np.vstack([baseData, row])
        i=i+1      
    
    columnNames=["MW","LogP","NumHDonors","NumHAcceptors"]   
    descriptors = pd.DataFrame(data=baseData,columns=columnNames)
    
    return descriptors


def pIC50(input):
    pIC50 = []

    for i in input['standard_value_norm']:
        molar = i*(10**-9) # Converts nM to M
        pIC50.append(-np.log10(molar))

    input['pIC50'] = pIC50
    x = input.drop('standard_value_norm', axis=1)
        
    return x

def norm_value(input):
    norm = []

    for i in input['standard_value']:
        if i > 100000000:
          i = 100000000
        norm.append(i)

    input['standard_value_norm'] = norm
    x = input.drop('standard_value', axis=1)
        
    return x


def exploratory_data(archivo_csv):
    df = pd.read_csv(archivo_csv)
    df_lipinski = lipinski(df.canonical_smiles)
    df_combined = pd.concat([df,df_lipinski], axis=1)
    df_norm = norm_value(df_combined)
    df_final = pIC50(df_norm)
    df_2class = df_final[df_final['class'] != 'intermediate']
    df_2class.to_csv('bioactivity_data_pi50.csv', index=False)

def exploratory_data1(archivo_csv):
    df = pd.read_csv(archivo_csv)

    df = pd.concat([df, lipinski(df.canonical_smiles)], axis=1)
    df = norm_value(df)
    df = pIC50(df)

    df[df['class'] != 'intermediate'].to_csv('bioactivity_data_pi50.csv', index=False)
