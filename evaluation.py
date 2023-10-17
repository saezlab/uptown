import pandas as pd
import os
import networkx as nx
import ptitprince as pt
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

class Eval:
    def __init__(self, dirpath, study_id):
        self.dirpath = dirpath
        self.G = None
        self.study_id = study_id
        self.target_dict = None
        self.tests_dict = self.get_test_dict()
        self.graphs, self.graphdata_df = self.parse_sif_files()
        self.distance_df = None
        self.filtered_df = None

    def get_datafiles(self):
        files_list = os.listdir(self.dirpath)
        sif_files = [os.path.join(self.dirpath, f) for f in files_list if (f.endswith('.sif') and (self.study_id in f))]
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
        iter_ids = []
        methods_list = []
        pagerank_thresholds = []

        file_paths = self.get_datafiles()
        
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
            iter_ids.append(parts[1])
            
            method_parts = parts[2:]
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
            'Iteration': iter_ids,
            'Methods': methods_list,
            'Pagerank Threshold': pagerank_thresholds
        })

        metadata_df.sort_values(by=['Study ID', 'Iteration', 'Methods', 'Pagerank Threshold'], inplace=True)
        
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

    def get_connected_targets(self, target_dict):
        """
        Given a network and a set of targets, compute how many targets are in the network (nodes)
        """
        num_targets = []
        for graph in self.graphs:
            iter_id = graph.split('__')[1]
            drug = iter_id.split('_')[1]
            targets = target_dict[drug].keys()
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
        iterations = set(self.graphdata_df['Iteration'].tolist())
        filtered_df = pd.DataFrame()

        for iterat in iterations:
            for method in methods:
                print(f"Processing {method} {iterat}...")
                sorted_df = self.graphdata_df[(self.graphdata_df['Methods'] == method) & (self.graphdata_df['Iteration'] == iterat)].sort_values(by=['Pagerank Threshold'])

                targets = self.target_dict[iterat].keys()
                num_targets = len(targets)
                # add column of missing targets
                sorted_df['Perc missing targets'] = (num_targets - sorted_df['Connected targets']) / num_targets * 100

                perc_missing_targets = sorted_df['Perc missing targets'].tolist()
                num_edges = sorted_df['Number of edges'].tolist()
                
                selected_threshold_index = self.find_elbow(num_edges, perc_missing_targets)
                # bind rows the selected threshold to the filtered df using pd.concat
                filtered_df = pd.concat([filtered_df, sorted_df.iloc[selected_threshold_index:selected_threshold_index+1, :]])
        
        return filtered_df
   



