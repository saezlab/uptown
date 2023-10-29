
import pandas as pd
from evaluation import Eval
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx

network_file = 'network_collectri.sif'

G = nx.read_weighted_edgelist(network_file, delimiter = '\t', create_using = nx.DiGraph)
for u, v, data in G.edges(data=True):
    weight = data['weight']
    data['sign'] = 1 if weight >= 0 else -1
    data['weight'] = abs(weight)

offtargets_df = pd.read_csv('panacea_offtargets.tsv', sep='\t')
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


G_eval = Eval(G, './results/', 'PANACEA')
G_eval.graphdata_df
G_eval.get_number_nodes()
G_eval.get_number_edges()
G_eval.get_connected_targets()
G_eval.threshold_filter()
G_eval.compute_centrality_metrics()
G_eval.get_returned_offtargets(offtarget_dict)

G_eval.graphdata_df.to_csv('PANACEA_graphdata_df.csv', index=False)
