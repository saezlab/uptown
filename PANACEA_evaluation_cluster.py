import pandas as pd
from evaluation import Eval
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import argparse
from path_calc import Solver

parser = argparse.ArgumentParser(description='PANACEA analysis')
parser.add_argument('-n', '--network', help='network file', required=True)
parser.add_argument('-s', '--sources', help='sources file', required=True) #panacea_sources.tsv
parser.add_argument('-t', '--tf', help='tf activity file', required=True) #tf_activity_results.tsv
parser.add_argument('-r', '--random', help='random_label', required=True)
parser.add_argument('-b', '--bio_context', help='combination cell_line-drug', required=True)
parser.add_argument('-d', '--dirpath', help='directory path', required=True)
parser.add_argument('-o', '--offtargets', help='offtargets file', required=True)
parser.add_argument('-i', '--dataset', help='dataset name', required=False, default='PANACEA')
parser.add_argument('-p', '--phospho', help='phospho file', required=False, default='')
parser.add_argument('-e', '--ec50_file', help='ec50 file', required=False, default='')

args = parser.parse_args()

network_file = args.network
tf_file = args.tf
random_label = args.random
bio_context = args.bio_context
dirpath = args.dirpath
offtargets_file = args.offtargets
dataset = args.dataset
phospho_file = args.phospho
ec50_file = args.ec50_file
sources_file = args.sources


G = nx.read_weighted_edgelist(network_file, delimiter = '\t', create_using = nx.DiGraph)
for u, v, data in G.edges(data=True):
    weight = data['weight']
    data['sign'] = 1 if weight >= 0 else -1
    data['weight'] = abs(weight)

offtargets_df = pd.read_csv(offtargets_file, sep='\t')
offtarget_dict = {}

# source
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

# Iterate over each row in the DataFrame
for index, row in offtargets_df.iterrows():
    # Get the drug and target from the row
    drug = row['cmpd']
    target = row['target']
    
    # Check if the drug is already a key in the dictionary
    if drug in offtarget_dict:
        # If it is, append the target to the list of targets for that drug
        offtarget_dict[drug].append(target)
    else:
        # If it isn’t, create a new list with the target and assign it to the drug key
        offtarget_dict[drug] = [target]

phospho_df = pd.read_csv(phospho_file, sep='\t')
phospho_dict = {}

ec50_df = pd.read_csv(ec50_file, sep='\t')
# filter ec50_df by biocontext (cellline_drug)
cell_line, drug = bio_context.split('_')
ec50_df = ec50_df[(ec50_df['Cell_line'] == cell_line) & (ec50_df['Drug'] == drug)]

for index, row in phospho_df.iterrows():
    # Get the drug and target from the row
    drug = row['drug']
    target = row['Gene']
    
    # Check if the drug is already a key in the dictionary
    if drug in phospho_dict:
        # If it is, append the target to the list of targets for that drug
        phospho_dict[drug].append(target)
    else:
        # If it isn’t, create a new list with the target and assign it to the drug key
        phospho_dict[drug] = [target]


G_eval = Eval(G, dirpath, dataset, bio_context, random_label)
# if there are no graphs in dictionary, finish execution
if not G_eval.graphs:
    exit()

cell_line, drug = bio_context.split('_')
try:
    G_eval.source_dict = source_dict[drug]
except KeyError:
    print('No sources for drug {}'.format(drug))
    exit()

G_eval.G = G_eval.reachability_filter(G_eval.G)

# filter the offtargets dictionary by the targets in the network
offtarget_dict = {k: v for k, v in offtarget_dict.items() if k in G_eval.G.nodes}

# filter the phospho df by the targets in the network
ec50_df = ec50_df[ec50_df['Gene'].isin(G_eval.G.nodes)]

G_eval.graphdata_df
G_eval.get_number_nodes()
G_eval.get_number_edges()
G_eval.get_connected_targets()
G_eval.threshold_filter()
G_eval.compute_centrality_metrics()
G_eval.get_returned_offtargets(offtarget_dict, 'offtargets')
G_eval.get_returned_offtargets(phospho_dict, 'phospho')
G_eval.compute_ec50_metrics(ec50_df)

G_eval.graphdata_df.to_csv(f'{dataset}_{random_label}_{bio_context}_graphdata_df.csv', index=False)
