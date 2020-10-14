import numpy as np
import pandas as pd
import os
from pathlib import Path


def parser(file_list, base_dir, delimiter, rank='species'):
    """
    Parses Kraken2 output files according to taxonomic rank of choice.

    Args:
        list containing paths of output files (created using --report-files in krakenuniq),
        base directory containing reports
        taxonomic rank of choice


    returns: dataframe with runid as row index and taxons as column index
    """
    # For testing
    # base_dir = Path.cwd() / 'datasets/kapusta_reports'
    # # file_list = os.listdir(base_dir)
    # # file_name = file_list[1]
    # rank = 'G'

    problem_files = []
    df_full = pd.DataFrame(columns=[rank])
    for file_name in file_list:
        try:
            df = pd.read_csv(base_dir / file_name, sep='\t', header=None)
            df.columns = ['%', 'cum_reads', 'reads', 'rank', 'taxID', 'taxName']
            runid = file_name.split(delimiter)[0]

            # Retrieve taxa rank of choice
            assert rank in df['rank'].unique()
            rank_df = df[df['rank'] == rank]
            rank_df = rank_df[['taxName', 'cum_reads']]

            # Strip whitespace
            rank_df.loc[:, 'taxName'] = rank_df.loc[:, 'taxName'].str.strip()

            # Drop unknown taxa
            rank_df = rank_df.loc[rank_df.taxName != 'uncultured', :]
            rank_df = rank_df.loc[rank_df.taxName != 'Unknown Family', :]
            rank_df = rank_df.loc[rank_df.taxName != 'Incertae Sedis', :]

            # Check for duplicate genera
            assert rank_df.duplicated().sum() == 0

            # Rename columns
            rank_df.columns = [rank, runid]
            rank_df.loc[:, runid] = rank_df.loc[:, runid].astype('int32')

            # # Set taxa as index
            # rank_df = rank_df.set_index(rank)

            # Outer join to get most taxa
            df_full = df_full.merge(rank_df, how='outer', on=rank)

        except AssertionError:
            problem_files.append(file_name)

        except ValueError:
            print(file_name)

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
reports = cwd / 'datasets/kraken2_reports'

# Kapusta
db = 'silva'
rank = 'G'

kapusta_list = os.listdir(reports / f'kapusta_reports_{db}')

# Parse reports
df = parser(kapusta_list, reports / f'kapusta_reports_{db}', delimiter='_', rank=rank)

# Add metadata
meta = pd.read_csv(datasets / 'kapusta_metadata_290620.tsv', sep='\t')
meta = meta[['SAMPLE_ID', 'GROUP']].copy()
meta = meta.rename({'SAMPLE_ID': 'run', 'GROUP': 'y'}, axis=1)
meta = meta.set_index('run')

df = meta.join(df, how='inner')
df = df.replace({'control': 'healthy', 'sample': 'septic'})

# Drop zero columns
df = df.loc[:, df.any(0)]   # Remove zero columns

# Clean column names to fit xgboost requirements
df.columns = df.columns.str.replace('>', '')
df.columns = df.columns.str.replace(']', '')
df.columns = df.columns.str.replace('[', '')

df.to_csv(datasets / f'kapusta_genus_raw_{db}.csv', index=False)
############################## Grumaz ##################################################################################
# Parse reports
grumaz_list = os.listdir(reports / 'grumaz_reports_maxi')

df2 = parser(grumaz_list, reports / 'grumaz_reports_maxi', delimiter='.', rank='G')

# Add metadata
meta = pd.read_csv(datasets / 'grumaz_pooled_metadata.csv', sep=',')
meta = meta[['run_accession', 'y']].copy()
meta = meta.rename({'run_accession': 'run'}, axis=1)
meta = meta.set_index('run')

df2 = meta.join(df2, how='inner')

# Drop zero columns
df2 = df2.loc[:, df2.any(0)]   # Remove zero columns

# Clean column names to fit xgboost requirements
df2.columns = df2.columns.str.replace('>', '')
df2.columns = df2.columns.str.replace(']', '')
df2.columns = df2.columns.str.replace('[', '')

# df2.to_csv(datasets / 'grumaz_genus_raw_maxi.csv', index=False)

#################################### Karius ############################################################################

# Parse reports
karius_list = os.listdir(reports / 'karius_reports_maxi')

df3 = parser(karius_list, reports / 'karius_reports_maxi', delimiter='.', rank='G')

# Add metadata
meta = pd.read_csv(datasets / 'karius_parsed_metadata.csv', sep=',')
meta = meta[['run', 'pathogen', 'y']].copy()
meta = meta.set_index('run')

df3 = meta.join(df3, how='inner')

# Drop zero columns
df3 = df3.loc[:, df3.any(0)]   # Remove zero columns

# Clean column names to fit xgboost requirements
df3.columns = df3.columns.str.replace('>', '')
df3.columns = df3.columns.str.replace(']', '')
df3.columns = df3.columns.str.replace('[', '')

df3.to_csv(datasets / 'karius_genus_raw_maxi.csv', index=False, header=True)

# Merge datasets
assert db == 'silva' and rank == 'G'
df.insert(value='kapusta', loc=0, column='dataset')
df2.insert(value='grumaz', loc=0, column='dataset')
df3.insert(value='karius', loc=0, column='dataset')

final_df = pd.concat([df, df2, df3], axis=0, join='inner', ignore_index=True)
final_df.to_csv(datasets / 'kapusta_grumaz_karius_genus_raw.csv', sep=',', index=False, header=True)

