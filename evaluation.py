import pandas as pd
import os
import networkx as nx
import ptitprince as pt
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import itertools
from path_calc import Solver    
# import umap
import matplotlib.pyplot as plt

class Eval:
    def __init__(self, G, dirpath, study_id, biocontext, random_label):
        self.dirpath = dirpath
        self.G = G
        self.study_id = study_id
        self.biocontext = biocontext
        self.random_label = random_label
        # self.tests_dict = self.get_test_dict()
        self.graphs, self.graphdata_df = self.parse_sif_files()
        self.target_dict = self.get_target_dict()
        self.distance_df = None
        __super__ = Solver
    
    def get_target_dict(self):
        self.target_dict = {}
        try:
            filename = f'./{self.dirpath}/{self.study_id}_{self.random_label}_{self.biocontext}_selectedtargets.csv'
            # convert df to dictionary and append to target dict
            df = pd.read_csv(filename, index_col=0)
            
            target_dict = df.to_dict()
            self.target_dict.update(target_dict)
            return self.target_dict
        except FileNotFoundError:
            return



    def get_datafiles(self):
        files_list = os.listdir(self.dirpath)
        sif_files = [os.path.join(self.dirpath, f) for f in files_list if (f.endswith('.sif') and (self.study_id in f) and (self.random_label in f) and (self.biocontext in f))]
        return sif_files
    
    def get_test_dict(self):
        filename = f'./{self.study_id}__testtargets.csv'
        try:
            df = pd.read_csv(filename, index_col=0)
            data_dict = df.to_dict()
            for key in data_dict:
                data_dict[key] = {k: v for k, v in data_dict[key].items() if pd.notna(v)}
            
            return data_dict
        
        except FileNotFoundError:
            
            return f"File {filename} not found."

    def parse_sif_files(self):
        """
        Parse a list of SIF files, return a dictionary of networkx graphs, and extract metadata with ordered methods.
        
        Parameters:
        - file_paths (list): List of paths to the SIF files.
        
        Returns:
        - dict: Dictionary where keys are filenames (without .sif extension) and values are networkx graphs.
        - DataFrame: Metadata dataframe with study ID, methods, pagerank threshold, and original graph ID columns.
        """
        graphs = {}
        
        # Lists to store the metadata
        graph_ids = []
        study_ids = []
        random_ids = []
        biocontext_ids = []
        methods_list = []
        pagerank_thresholds = []

        file_paths = self.get_datafiles()
        # if it is emmpty, finish execution
        if not file_paths:
            return None, None
        
        for file_path in file_paths:
            # Extract the filename without the .sif extension
            graph_id = file_path.split('/')[-1].replace('.sif', '')
            
            # Store the original graph ID
            graph_ids.append(graph_id)
            
            # Initialize a directed graph
            G = nx.DiGraph()
            
            # Parse the SIF file, skipping the header
            with open(file_path, 'r') as f:
                lines = f.readlines()[1:]  # skip the header
                for line in lines:
                    parts = line.strip().split('\t')
                    if len(parts) == 3:
                        source, interaction, target = parts
                        interaction = 1 if interaction == 'P' else (-1 if interaction == 'N' else interaction)
                        G.add_edge(source, target, interaction=interaction)
            
            # Store the graph in the dictionary
            graphs[graph_id] = G
            
            # Metadata extraction
            parts = graph_id.split('__')
            study_ids.append(parts[0])
            random_ids.append(parts[1])
            biocontext_ids.append(parts[2])
            
            method_parts = parts[3:]
            pagerank_value = None
            for idx, part in enumerate(method_parts):
                if "pagerank" in part:
                    pagerank_value = float(part.split('_')[-1])  # Extract float value after "pagerank_"
                    method_parts[idx] = "pagerank"  # Replace with "pagerank" string
            
            pagerank_thresholds.append(pagerank_value)
            methods_list.append(', '.join(method_parts))
        
        # Create the metadata dataframe
        metadata_df = pd.DataFrame({
            'Graph ID': graph_ids,
            'Study ID': study_ids,
            'Random': random_ids,
            'Biocontext': biocontext_ids,
            'Methods': methods_list,
            'Pagerank Threshold': pagerank_thresholds
        })

        metadata_df.sort_values(by=['Study ID', 'Random', 'Biocontext', 'Methods', 'Pagerank Threshold'], inplace=True)
        
        return graphs, metadata_df
    

    # Distance eval
    def distance_calc(self):
        distances = []
        test_nodes = []
        network_nodes = []
        graph_ids = []
        run_ids = []
        undirG = self.G.reverse()
        for graph in self.graphs:
            run = graph.strip().split('__')[1]
            nodes_subnetwork = [node for node in self.graphs[graph].nodes]
            for node in nodes_subnetwork:
                for test_node in self.tests_dict[run]:
                    try:
                        distance = nx.shortest_path_length(undirG, source=test_node, target=node)
                    except nx.NetworkXNoPath:
                        distance = np.Inf
                    distances.append(distance)
                    test_nodes.append(test_node)
                    network_nodes.append(node)
                    graph_ids.append(graph)
                    run_ids.append(run)
        
        self.distance_df = pd.DataFrame({
            'Graph ID': graph_ids,
            'Run ID': run_ids,
            'Study ID': self.study_id,
            'Test nodes': test_nodes,
            'Network nodes': network_nodes,
            'Distances': distances
        })

        self.plot_distance(self.distance_df)
        
        return self.distance_df
    
    def plot_distance(self, distances_df):
        distances_df["Main Method"] = distances_df["Graph ID"].apply(lambda x: x.split("__")[2].split('_pagerank')[0])
        plt.figure(figsize=(15, 10))
        
        # Creating a color palette for pagerank thresholds
        unique_thresholds = distances_df["Graph ID"].apply(lambda x: x.split('_')[-1]).unique()
        color_palette = dict(zip(unique_thresholds, sns.color_palette("rainbow", len(unique_thresholds))))
        
        # Boxplot
        sns.boxplot(data=distances_df, x="Main Method", y="Distances", hue="Graph ID", palette=color_palette, dodge=True)
        
        # Adjusting the plot aesthetics
        plt.xticks(rotation=45)
        plt.legend(title="Pagerank Threshold", bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.tight_layout()
        
        plt.show()

    # Comparisons
    def get_number_nodes(self):
        """
        For a given dictionary of networks, computes the number of nodes and appends it to the self.graphdata_df dataframe.
        """
        num_nodes = []
        for graph in self.graphs:
            num_nodes.append(len(self.graphs[graph].nodes))
        nodes_info = pd.DataFrame({'Graph ID': list(self.graphs.keys()), 'Number of nodes': num_nodes})
        self.graphdata_df = pd.merge(self.graphdata_df, nodes_info, on="Graph ID")
    
    def get_number_edges(self):
        num_edges = []
        for graph in self.graphs:
            num_edges.append(len(self.graphs[graph].edges))
        edges_info = pd.DataFrame({'Graph ID': list(self.graphs.keys()), 'Number of edges': num_edges})
        self.graphdata_df = pd.merge(self.graphdata_df, edges_info, on="Graph ID")

    def get_connected_targets(self):
        """
        Given a network and a set of targets, compute how many targets are in the network (nodes)
        """
        num_targets = []
        for graph in self.graphs:
            parts = graph.split('__')
            iter_id = '__'.join(parts[0:4])
            targets = self.target_dict[iter_id].keys()
            num_targets.append(len(set(targets) & set(self.graphs[graph].nodes)))
        targets_info = pd.DataFrame({'Graph ID': list(self.graphs.keys()), 'Connected targets': num_targets})
        self.graphdata_df = pd.merge(self.graphdata_df, targets_info, on="Graph ID")

    def compute_degree_distribution(self):
        """
        Compute the degree distribution for each graph in the self.graphs dictionary.
        """
        degree_distributions = []
        for graph in self.graphs:
            degree_distributions.append(self.graphs[graph].degree)
        degree_info = pd.DataFrame({'Graph ID': list(self.graphs.keys()), 'Degree distribution': degree_distributions})
        self.graphdata_df = pd.merge(self.graphdata_df, degree_info, on="Graph ID")

    def find_elbow(self, x, y):
        """
        Find the elbow of a curve using x and y coordinates.

        Args:
            x (list): List of x coordinates.
            y (list): List of y coordinates.

        Returns:
            int: Index of the elbow point.
        """
        # Coordinates of the first point
        p1 = np.array([x[0], y[0]])
        
        # Coordinates of the last point
        p2 = np.array([x[-1], y[-1]])
        
        
        # Find the elbow by computing distance from each point to the line
        distances = []
        for i in range(len(x)):
            point = np.array([x[i], y[i]])
            distance = np.linalg.norm(np.cross(p2-p1, p1-point))/np.linalg.norm(p2-p1)
            distances.append(distance)
            
        # Return the index of the point with max distance which is the elbow
        return np.argmax(distances)
    
    def threshold_filter(self):
        perc_missing_targets = []
        num_edges = []

        methods = set(self.graphdata_df['Methods'].tolist())
        random_ids = set(self.graphdata_df['Random'].tolist())
        biocontexts = set(self.graphdata_df['Biocontext'].tolist())
        filtered_df = pd.DataFrame()
        
        for random_id in random_ids:
            for biocontext in biocontexts:
                for method in methods:
                    drug = biocontext.split('_')[1]
                    cell_line = biocontext.split('_')[0]
                    # print(f"Processing {method} {biocontext} {random_id}...")
                    sorted_df = self.graphdata_df[(self.graphdata_df['Methods'] == method) & (self.graphdata_df['Biocontext'] == biocontext) & (self.graphdata_df['Random'] == random_id)].sort_values(by=['Pagerank Threshold'])
                    parts = set(sorted_df['Graph ID'].tolist())
                    iter_id = '__'.join(parts.pop().split('__')[0:4])

                    targets = self.target_dict[iter_id].keys()
                    num_targets = len(targets)
                    # print(num_targets)
                    
                    # add column of missing targets
                    sorted_df['Perc missing targets'] = (num_targets - sorted_df['Connected targets']) / num_targets * 100

                    perc_missing_targets = sorted_df['Perc missing targets'].tolist()
                    num_edges = sorted_df['Number of edges'].tolist()
                    
                    selected_threshold_index = self.find_elbow(num_edges, perc_missing_targets)
                    # bind rows the selected threshold to the filtered df using pd.concat
                    threshold_100 = sorted_df[sorted_df['Pagerank Threshold'] == 100]
                    filtered_df = pd.concat([filtered_df, sorted_df.iloc[selected_threshold_index:selected_threshold_index+1, :], threshold_100])
            
        filtered_df = filtered_df[~((filtered_df['Pagerank Threshold'] == 100) & (filtered_df['Methods'].str.contains('reachability, pagerank, reachability, allpaths')))]
        
        filtered_df.loc[filtered_df['Pagerank Threshold'] == 100, 'Methods'] = filtered_df['Methods'].apply(lambda x: x.replace('pagerank, ', ''))
        filtered_df['cell_line'] = filtered_df['Biocontext'].apply(lambda x: x.split('_')[0]).tolist()
        filtered_df['drug'] = filtered_df['Biocontext'].apply(lambda x: x.split('_')[1]).tolist()

        
        self.graphdata_df = filtered_df
    
    def get_returned_offtargets(self, offtarget_dict):
        offtarget_results = []

        # Iterating over each key in the network dictionary
        for key, graph in self.graphs.items():
            # check if the graph is in the graphdata_df
            if key not in self.graphdata_df['Graph ID'].tolist():
                continue
            # Iterating over each drug and its targets in the drug_target_dict
            number_nodes = len(graph.nodes)
            number_edges = len(graph.edges)
            for drug, targets in offtarget_dict.items():
                # Checking if the drug name is present in the key
                if drug.lower() in key.lower():
                    # Counting how many targets are present in the nodes of the network
                    target_count = sum(1 for target in targets if target in graph.nodes)
                    all_offtargets = len(targets)
                    perc_offtargets = (target_count / len(targets)) * 100
                    perc_offtargets_nodes = (target_count / number_nodes) * 100
                    perc_offtargets_edges = (target_count / number_edges) * 100
                    # Appending the results to the list
                    offtarget_results.append({'Graph ID': key, 'offtarget_count': target_count, 'all_offtargets': all_offtargets, 'perc_offtarget': perc_offtargets, 'perc_offtarget_nodes': perc_offtargets_nodes, 'perc_offtarget_edges': perc_offtargets_edges})
            # if the drug is not present in the offtarget_dict, append NA to the offtarget_results
            if not any(drug.lower() in key.lower() for drug in offtarget_dict.keys()):
                offtarget_results.append({'Graph ID': key, 'offtarget_count': np.nan, 'all_offtargets': np.nan, 'perc_offtarget': np.nan, 'perc_offtarget_nodes': np.nan, 'perc_offtarget_edges': np.nan})
        # Creating a DataFrame from the results list
        offtarget_result_df = pd.DataFrame(offtarget_results)
        self.graphdata_df = pd.merge(self.graphdata_df, offtarget_result_df, on="Graph ID")
        
    def compute_overlap(self, iter):
        subset_df = self.graphdata_df[self.graphdata_df['Iteration'] == iter]
        print(subset_df)
    
    def degree_centrality(self):
        degree_centrality = []

        for graph in self.graphs:
            degree_results = nx.degree_centrality(self.graphs[graph])
            degree_centrality.append({'Graph ID': graph, 
                                      'Degree centrality': degree_results, 
                                      'Mean degree centrality': np.mean(list(degree_results.values()))})




        degree_centrality_df = pd.DataFrame(degree_centrality)
        self.graphdata_df = pd.merge(self.graphdata_df, degree_centrality_df, on="Graph ID")


    def closeness_centrality(self):
        closeness_centrality = []

        for graph in self.graphs:
            closeness_results = nx.closeness_centrality(self.graphs[graph])
            closeness_centrality.append({'Graph ID': graph, 
                                         'Closeness centrality': closeness_results, 
                                         'Mean closeness centrality': np.mean(list(closeness_results.values()))})
        closeness_centrality_df = pd.DataFrame(closeness_centrality)
        self.graphdata_df = pd.merge(self.graphdata_df, closeness_centrality_df, on="Graph ID")
    
    def betweenness_centrality(self):
        betweenness_centrality = []

        for graph in self.graphs:
            betweenness_results = nx.betweenness_centrality(self.graphs[graph])
            betweenness_centrality.append({'Graph ID': graph, 
                                           'Betweenness centrality': betweenness_results,
                                           'Mean betweenness centrality': np.mean(list(betweenness_results.values()))})
        betweenness_centrality_df = pd.DataFrame(betweenness_centrality)
        self.graphdata_df = pd.merge(self.graphdata_df, betweenness_centrality_df, on="Graph ID")
    
    def compute_centrality_metrics(self):
        print('Computing degree centrality')
        self.degree_centrality()
        print('Computing closeness centrality')
        self.closeness_centrality()
        print('Computing betweenness centrality')
        self.betweenness_centrality()

    def create_edge_matrix(self):
        """
        Creates a binary matrix representing the presence or absence of edges in each graph.
        
        Parameters:
        - graph_dict: dictionary
            Keys are graph IDs and values are NetworkX graphs.
        
        Returns:
        - pandas DataFrame
            Rows represent edges, columns represent graphs, values are binary.
        """
        # Extract all unique edges from all graphs
        all_edges = set()
        graph_ids = self.graphdata_df['Graph ID'].tolist()
        #get only graphs present in the graph ids list
        graph_dict = {k: v for k, v in self.graphs.items() if k in graph_ids}
        for graph in graph_dict.values():
            all_edges.update(graph.edges)
        
        # Sort edges for consistent ordering
        all_edges = sorted(all_edges)
        
        # Create the matrix
        matrix_data = []
        for edge in all_edges:
            row = []
            for graph_id, graph in graph_dict.items():
                row.append(1 if graph.has_edge(*edge) else 0)
            matrix_data.append(row)
        
        # Create DataFrame
        columns = list(graph_dict.keys())
        edge_matrix = pd.DataFrame(matrix_data, columns=columns, index=all_edges)
        
        return edge_matrix
    
    def plot_umap(self, features):
        edge_matrix = self.create_edge_matrix()
        edge_metric_columns = edge_matrix.columns.tolist()
        filtered_rundf = self.graphdata_df
        # set Graph ID as index
        filtered_rundf = filtered_rundf.set_index('Graph ID')
        # sort the filtered_rundf by the order of the edge metric column ids and the Graph ID column from the filtered_rundf
        filtered_rundf = filtered_rundf.reindex(index=edge_metric_columns)
        # Applying UMAP to the edge matrix
        reducer = umap.UMAP()
        embedding = reducer.fit(edge_matrix.transpose())

        for feature in features:
            umap.plot.points(embedding, labels=filtered_rundf[feature])
    
    def plot_raincloud_plots(self, feature):
        plt.figure(figsize=(15, 10))
        ax = pt.RainCloud(data = self.graphdata_df, x = 'Methods', y = feature, palette = "Set2",
                            orient = 'h', width_box= .1, width_viol=1)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45,  ha='right')

        plt.show()
