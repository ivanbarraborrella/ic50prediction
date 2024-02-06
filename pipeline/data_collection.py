import pandas as pd
from chembl_webresource_client.new_client import new_client

def data_collection(protein):
    target = new_client.target
    target_query = target.search(protein)
    targets = pd.DataFrame.from_dict(target_query)
    selected_target = targets.target_chembl_id[0]

    activity = new_client.activity
    res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
    df = pd.DataFrame.from_dict(res)

    df = df[df.standard_value.notna() & df.canonical_smiles.notna()]
    df = df.drop_duplicates(['canonical_smiles'])

    selection = ['molecule_chembl_id','canonical_smiles','standard_value']
    df = df[selection]

    df['class'] = df['standard_value'].apply(lambda x: 'active' if float(x) <= 1000 else ('inactive' if float(x) >= 10000 else 'intermediate'))

    df.to_csv('bioactivity_data_curated.csv', index=False)