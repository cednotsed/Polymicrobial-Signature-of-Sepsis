from pathlib import Path
import numpy as np
import pandas as pd
import os
cwd = Path.cwd()
datasets = cwd / 'datasets'
metadata = datasets / 'metadata'

exp = pd.read_csv(metadata / 'experiment_name_karius.csv', dtype='str')
exp = exp[['Experiment Accession', 'Experiment Title']]
exp.columns = ['exp_id', 'exp_name']
exp_types = exp['exp_name'].unique()

pos = exp[exp['exp_name'] == exp_types[0]]
neg = exp[exp['exp_name'] == exp_types[1]]
combined = pd.concat([pos, neg], axis=0)

runtable = pd.read_csv(metadata / 'SraRunTable.txt', dtype='str')
runtable = runtable[['Run', 'Experiment', 'subject_id', 'Sample Name']]
runtable.columns = ['run', 'exp_id', 'subject_id', 'sample_name']
final = pd.merge(runtable, combined, on='exp_id', how='inner')

definite = pd.read_csv(metadata / 'definite_id_karius.csv', dtype='str')
definite = definite[['Subject ID', 'Karius Result', 'Confimatory Test']]
definite.columns = ['subject_id', 'pathogen', 'test']

septic = pd.merge(definite, final, how='inner', on='subject_id')
septic = septic.dropna(axis=0, subset=['test'])

healthy = final[final['exp_name'] == exp_types[1]]
healthy.insert(loc=0, column='pathogen', value='none')
healthy.insert(loc=0, column='test', value='none')
df = pd.concat([healthy, septic], sort=True, axis=0)
df = df.replace({exp_types[0]: 'septic', exp_types[1]: 'healthy'})
df.to_csv(metadata / 'parsed_metadata.csv', index=False)
