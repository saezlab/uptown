import pandas as pd
from solver import Solver
import networkx as nx
import argparse

# Parameters

# Helper functions
def custom_sorting_dict(df, bio_context):
    target_dict = {}
    df['abs_values'] = df[bio_context].apply(abs)

    # Sorting the DataFrame based on the temporary column
    df.sort_values(by='abs_values', inplace=True, ascending=False)

    # Dropping the temporary column
    df.drop(columns=['abs_values'], inplace=True)

    bio_context_dict = df[bio_context].to_dict()

    target_dict[bio_context] = bio_context_dict
    
    return target_dict

# MENU
parser = argparse.ArgumentParser(description='PANACEA analysis')
parser.add_argument('-n', '--network', help='network file', required=True)
parser.add_argument('-s', '--sources', help='sources file', required=True) #panacea_sources.tsv
parser.add_argument('-t', '--tf', help='tf activity file', required=True) #tf_activity_results.tsv
parser.add_argument('-r', '--random', help='random_label', required=True)
parser.add_argument('-b', '--bio_context', help='combination cell_line-drug', required=True)
parser.add_argument('-tf', '--tfs', help='number of tfs to use', required=False, default=25)
parser.add_argument('-a', '--allpaths', help='allpaths cutoff', required=False, default=6)
parser.add_argument('-d', '--dataset', help='dataset name', required=False, default='PANACEA')

args = parser.parse_args()

network_file = args.network
sources_file = args.sources
tf_file = args.tf
random_label = args.random
number_tfs = int(args.tfs)
bio_context = args.bio_context
cutoff = int(args.allpaths)
dataset = args.dataset


# Network import
G = nx.read_weighted_edgelist(network_file, delimiter = '\t', create_using = nx.DiGraph)
for u, v, data in G.edges(data=True):
    weight = data['weight']
    data['sign'] = 1 if weight >= 0 else -1
    data['weight'] = abs(weight)

## Sources (KI targets) formatting
source_df = pd.read_csv(sources_file, sep='\t')

nodes_network = [f for f in G.nodes]
filtered_source_df = source_df[source_df.target.isin(nodes_network)]

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
tf_activities = pd.read_csv(tf_file, sep='\t', index_col=0).transpose()

tf_dict = custom_sorting_dict(tf_activities, bio_context)

# Run the solver
selected_targets_dict = {} 

# Generate a random string of 5 characters
# not equal to real:
if random_label != 'real':
    random_part = random_label
    random_bool = True
elif random_label == 'real':
    random_part = 'real'
    random_bool = False

cell_line, drug = bio_context.split('_')
        
print('Solving for cell line and drug {}'.format(bio_context))    
G_solver = Solver(G, dataset)
G_solver.random = random_bool
try:
    G_solver.source_dict = source_dict[drug]
except KeyError:
    print('No sources for drug {}'.format(drug))
    exit()
G_solver.network_batchrun(iter = bio_context, 
                            random_ident = random_part, 
                            tf_dict = tf_dict, 
                            number_tfs=number_tfs,
                            cutoff = cutoff)
selected_targets_dict.update(G_solver.selected_targets_dict)

targets_df = pd.DataFrame.from_dict(selected_targets_dict)
targets_df.to_csv(f'{dataset}_{random_label}_{bio_context}_selectedtargets.csv')
    


