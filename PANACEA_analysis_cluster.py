import pandas as pd
from path_calc import Solver
import networkx as nx
import string
import random

# Parameters
random_iterations = 1000
number_tfs = 25

# Helper functions
def custom_sorting_dict(df, bio_contexts):
    target_dict = {}
    for bio_context in bio_contexts:
        df['abs_values'] = df[bio_context].apply(abs)

        # Sorting the DataFrame based on the temporary column
        df.sort_values(by='abs_values', inplace=True, ascending=False)

        # Dropping the temporary column
        df.drop(columns=['abs_values'], inplace=True)

        bio_context_dict = df[bio_context].to_dict()

        target_dict[bio_context] = bio_context_dict
    
    return target_dict

# Network import
G = nx.read_weighted_edgelist('network_collectri.sif', delimiter = '\t', create_using = nx.DiGraph)
for u, v, data in G.edges(data=True):
    weight = data['weight']
    data['sign'] = 1 if weight >= 0 else -1
    data['weight'] = abs(weight)

## Sources (KI targets) formatting
source_df = pd.read_csv('panacea_sources.tsv', sep='\t')

nodes_network = [f for f in G.nodes]
filtered_source_df = source_df[source_df.target.isin(nodes_network)]
filtered_source_df

# create a dictionary per treatment, with the targets as keys and the sign as value
source_dict = {}
for i in range(len(filtered_source_df)):
    treatment = filtered_source_df.iloc[i, 0]
    target = filtered_source_df.iloc[i, 1]
    sign = filtered_source_df.iloc[i, 2]
    if treatment not in source_dict:
        source_dict[treatment] = {}
    source_dict[treatment][target] = float(sign)

## TF formatting
tf_activities = pd.read_csv('tf_activity_results.tsv', sep='\t', index_col=0).transpose()
bio_contexts = tf_activities.columns.tolist()
tf_dict = custom_sorting_dict(tf_activities, bio_contexts)

# null model
random_list = [False] * 1 + [True] * random_iterations

# Run the solver
selected_targets_dict = {} 
for random_label in random_list:
    # Generate a random string of 5 characters
    if random_label:
        characters = string.ascii_letters + string.digits
        random_part = ''.join(random.choice(characters) for i in range(8))
    elif not random_label:
        random_part = 'real'

    for bio_context in bio_contexts:
        cell_line, drug = bio_context.split('_')
             
        print('Solving for cell line and drug {}'.format(bio_context))    
        G_solver = Solver(G, 'PANACEA')
        G_solver.random = random_label
        try:
            G_solver.source_dict = source_dict[drug]
        except KeyError:
            print('No sources for drug {}'.format(drug))
            continue
        G_solver.network_batchrun(iter = bio_context, 
                                  random_ident = random_part, 
                                  tf_dict = tf_dict, 
                                  number_tfs=number_tfs,
                                  cutoff = 6)
        selected_targets_dict.update(G_solver.selected_targets_dict)

targets_df = pd.DataFrame.from_dict(selected_targets_dict)
targets_df.to_csv('./results/PANACEA_selected_targets.csv')
    


