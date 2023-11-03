
import pandas as pd
from evaluation import Eval
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import argparse

parser = argparse.ArgumentParser(description='PANACEA analysis')
parser.add_argument('-n', '--network', help='network file', required=True)
parser.add_argument('-t', '--tf', help='tf activity file', required=True) #tf_activity_results.tsv
parser.add_argument('-r', '--random', help='random_label', required=True)
parser.add_argument('-b', '--bio_context', help='combination cell_line-drug', required=True)
parser.add_argument('-d', '--dirpath', help='directory path', required=True)
parser.add_argument('-o', '--offtargets', help='offtargets file', required=True)

args = parser.parse_args()

network_file = args.network
tf_file = args.tf
random_label = args.random
bio_context = args.bio_context
dirpath = args.dirpath
offtargets_file = args.offtargets



G = nx.read_weighted_edgelist(network_file, delimiter = '\t', create_using = nx.DiGraph)
for u, v, data in G.edges(data=True):
    weight = data['weight']
    data['sign'] = 1 if weight >= 0 else -1
    data['weight'] = abs(weight)

offtargets_df = pd.read_csv(offtargets_file, sep='\t')
offtarget_dict = {}

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


G_eval = Eval(G, dirpath, 'PANACEA', bio_context, random_label)
# if there are no graphs in dictionary, finish execution
if not G_eval.graphs:
    exit()
G_eval.graphdata_df
G_eval.get_number_nodes()
G_eval.get_number_edges()
G_eval.get_connected_targets()
G_eval.threshold_filter()
G_eval.compute_centrality_metrics()
G_eval.get_returned_offtargets(offtarget_dict)

G_eval.graphdata_df.to_csv(f'PANACEA_{random_label}_{bio_context}_graphdata_df.csv', index=False)
