import numpy as np
import pandas as pd
import os
from pathlib import Path


def parser(file_list, base_dir, rank='species'):
    """
    Parses krakenuniq output files according to taxonomic rank of choice.

    Args:
        list containing paths of output files (created using --report-files in krakenuniq),
        base directory containing reports
        taxonomic rank of choice


    returns: dataframe with runid as row index and taxons as column index
    """
    problem_files = []
    df_full = pd.DataFrame(columns=[rank])
    for file_name in file_list:
        try:
            # Column names ['%', 'reads', 'taxReads', 'kmers', 'dup', 'cov', 'taxID', 'rank', 'taxName']
            df = pd.read_csv(os.path.join(base_dir, file_name), sep='\t')
            runid = file_name.split(".")[0]

            # Retrieve taxa rank of choice
            assert rank in df['rank'].unique()
            rank_df = df[df['rank'] == rank]
            rank_df = rank_df[['taxName', 'kmers']]

            # # Normalise no. of kmers per rank
            # rank_df['kmers'] /= rank_df['kmers'].sum()

            # Rename columns
            rank_df.columns = [rank, runid]

            # Outer join to get most taxa
            df_full = df_full.merge(rank_df, how='outer', on=rank)

        except:
            problem_files.append(file_name)

    print("These files were in an invalid format:")
    print(*problem_files, sep='\n')

    df_full = df_full.set_index(rank)
    df_full = df_full.T

    # Replace Nans with zeros
    df_full = df_full.replace({np.nan: 0})

    return df_full


# Get Paths
cwd = Path.cwd()
datasets = cwd / 'datasets'
metadata = datasets / 'metadata'
karius_list = os.listdir(datasets / 'karius_reports')

# Parse Karius reports
karius_genus = parser(karius_list, datasets / 'karius_reports', rank='genus')
karius_genus = karius_genus.reset_index(drop=False)
karius_genus = karius_genus.rename({'index': 'run'}, axis=1)

karius_species = parser(karius_list, datasets / 'karius_reports', rank='species')
karius_species = karius_species.reset_index(drop=False)
karius_species = karius_species.rename({'index': 'run'}, axis=1)

# Add metadata
meta_df = pd.read_csv(metadata / 'karius_parsed_metadata.csv')
meta_df = meta_df[['pathogen', 'run', 'exp_name']]
meta_df = meta_df.drop_duplicates(subset=['run'])
meta_df.to_csv(datasets / 'sample_confirmed_list.csv', header=True, index=False)
print(set(meta_df['run']) - set(karius_genus['run']))
karius_genus = pd.merge(meta_df, karius_genus, how='inner', on='run')
karius_species = pd.merge(meta_df, karius_species, how='inner', on='run')

karius_genus = karius_genus.rename({'exp_name': 'y'}, axis=1)
karius_species = karius_species.rename({'exp_name': 'y'}, axis=1)

# Drop zero columns and clean column names
karius_genus = karius_genus.loc[:, karius_genus.any(0)]     # Remove zero columns
karius_genus.columns = karius_genus.columns.str.strip()

karius_species = karius_species.loc[:, karius_species.any(0)]
karius_species.columns = karius_species.columns.str.strip()

# Manipulating column names to fit xgboost requirements
karius_species.columns = karius_species.columns.str.replace('>', '')
karius_species.columns = karius_species.columns.str.replace(']', '')
karius_species.columns = karius_species.columns.str.replace('[', '')

# Get pathogen dataset
pathogens = pd.read_csv(metadata / 'pathogen_list.csv', encoding='ISO-8859–1')
pathogens = pathogens[pathogens['HostGroup'] == 'Human']
pathogenic_genera = set(pathogens['Genus'])
genera = set(karius_genus.columns)
intersecting_index = ['pathogen'] + ['y'] + list(pathogenic_genera.intersection(genera))
karius_genus.loc[:, intersecting_index]
karius_pathogens = karius_genus.loc[:, intersecting_index]
karius_genus = karius_genus.set_index(keys='run', drop=False)
karius_pathogens = karius_pathogens.set_index(keys=karius_genus.index)

import pickle
with open(datasets / 'karius_species_raw.pickle', 'wb') as f:
    pickle.dump((karius_species.iloc[:, 3:], karius_species['y']), f)

with open(datasets / 'karius_genus_raw.pickle', 'wb') as f:
    pickle.dump((karius_genus.iloc[:, 3:], karius_genus['y']), f)

with open(datasets / 'karius_genus_pathogens.pickle', 'wb') as f:
    pickle.dump((karius_pathogens.iloc[:, 2:], karius_pathogens['y']), f)

karius_genus.to_csv(datasets / 'karius_genus_raw.csv', index=False)
karius_species.to_csv(datasets / 'karius_species_raw.csv', index=False)
karius_pathogens.to_csv(datasets / 'karius_genus_pathogens.csv', index=False)
