import pandas as pd
import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
import os
from upsetplot import plot
import itertools
import time
import ptitprince as pt
import statannot as st
from pathos.multiprocessing import ProcessingPool as Pool


class Solver:
    """
    Solver for finding signed paths in a graph.
    """

    def __init__(self, G, study_id = None):
        """
        Initialize the signed path solver with a given graph.

        Args:
            G (networkx graph): A networkx graph object.
            study_id (Any, optional): An identifier for the study. Defaults to None.
        """
        self.G = G
        self.study_id = study_id
        self.label = 'pagerank'
        self.source_dict = {}
        self.target_dict = {}
        self.threshold = None
        self.subG = None
        self.is_reversed = False
        self.all_paths_res = []
        self.connected_all_path_targets = {}
        self.shortest_paths_res = []
        self.pagerank_values = {}
        self.connected_targets = {}
        self.sc_paths_res = []
        self.connected_sc_targets = {}
        self.runinfo_df = pd.DataFrame(columns=['label', 'num_nodes', 'num_edges', 'degrees', 'elapsed_time'])



    def get_subnetwork(self, paths, sign_consistent=True):
        """
        Creates a subnetwork from a list of paths. Also checks for sign consistency.

        Args:
            paths (list): A list of lists containing paths.
            sign_consistent (bool, optional): Determines if the network should check for sign consistency. Defaults to True.
        
        Returns:
            Tuple: Contains the subgraph, connected targets, and the resulting paths.
        """
        V = nx.DiGraph()
        connected_targets = {}
        paths_res = []
        for path in paths:
            source = path[0]
            product_sign = 1
            target = path[-1]
            source_sign = self.source_dict[source]
            target_sign = self.target_dict[target]

            if source not in connected_targets:
                connected_targets[source] = []

            for i in range(len(path) - 1):
                edge_sign = self.subG.get_edge_data(path[i], path[i + 1])['sign']
                
                product_sign *= edge_sign

            if (sign_consistent and source_sign * product_sign == target_sign) or not sign_consistent:
                paths_res.append(path)
                connected_targets[source].append(target) if target not in connected_targets[source] else None
                for i in range(len(path) - 1):
                    edge_data = self.subG.get_edge_data(path[i], path[i + 1])
                    V.add_edge(path[i], path[i + 1], **edge_data)
        
        return V, connected_targets, paths_res



    def pagerank_solver(self, alpha=0.85, max_iter=100, tol=1.0e-6, nstart=None, weight='weight', personalize_for="source"):
        """
        Compute the PageRank values for nodes in the graph.

        Args:
            alpha (float): Damping factor for the PageRank algorithm.
            max_iter (int): Maximum number of iterations.
            tol (float): Tolerance to determine convergence.
            nstart (dict): Starting value of PageRank iteration for all nodes.
            weight (str): Edge data key to use as weight.
            personalize_for (str): Personalize the PageRank by setting initial probabilities for either sources or targets.
        """
    
        if personalize_for == "source":
            personalized_prob = {n: 1/len(self.source_dict) for n in self.source_dict}
        elif personalize_for == "target":
            personalized_prob = {n: 1/len(self.target_dict) for n in self.target_dict}
        else:
            raise ValueError("personalize_for should be either 'source' or 'target'")
        
        
        if personalize_for == "target" and not self.is_reversed:
            self.subG = self.reverse_graph(self.subG)
        elif personalize_for == "source" and self.is_reversed:
            self.subG = self.reverse_graph(self.subG)

        personalized_prob.update({n: 0 for n in self.subG.nodes() if n not in personalized_prob})

        pagerank = nx.pagerank(self.subG, alpha=alpha, max_iter=max_iter, personalization=personalized_prob, tol=tol, nstart=nstart, weight=weight, dangling=personalized_prob)
        
        self.pagerank_values[personalize_for] = pagerank
        
        for node, pr_value in pagerank.items():
            attribute_name = 'pagerank_from_targets' if personalize_for == "target" else 'pagerank_from_sources'
            self.subG.nodes[node][attribute_name] = pr_value
        
        if personalize_for == "target" and self.is_reversed:
            self.reverse_graph(self.subG)



    def compute_overlap(self):
        """
        Compute the overlap of nodes that exceed the PageRank threshold from sources and targets.

        Returns:
            tuple: Contains nodes above threshold from sources, nodes above threshold from targets, and overlapping nodes.
        """
        nodes_above_threshold_from_sources = {node for node, data in self.subG.nodes(data=True) if data.get('pagerank_from_sources') > self.threshold}
        nodes_above_threshold_from_targets = {node for node, data in self.subG.nodes(data=True) if data.get('pagerank_from_targets') > self.threshold}
        overlap = nodes_above_threshold_from_sources.intersection(nodes_above_threshold_from_targets)
        nodes_to_include = nodes_above_threshold_from_sources.union(nodes_above_threshold_from_targets)
        self.subG = self.subG.subgraph(nodes_to_include)
        self.label = f'{self.label}__pagerank_{self.threshold}'
        return nodes_above_threshold_from_sources, nodes_above_threshold_from_targets, overlap



    def shortest_paths(self, verbose = False):
        """
        Calculate the shortest paths between sources and targets.

        Args:
            verbose (bool): If True, print warnings when no path is found to a given target.

        Returns:
            list: A list containing the shortest paths.
        """
        self.label = f'{self.label}__shortest'
        start_time = time.time()
        self.shortest_paths_res = []
        self.connected_targets = {}
        sources = list(self.source_dict.keys())
        targets = list(self.target_dict.keys())
        for source_node in sources:
            if source_node not in self.connected_targets:
                self.connected_targets[source_node] = []

            for target_node in targets:
                try:
                    self.shortest_paths_res.extend([p for p in nx.all_shortest_paths(self.subG, source=source_node, target=target_node, weight='weight')])
                    self.connected_targets[source_node].append(target_node)
                except nx.NetworkXNoPath as e:
                    if verbose:
                        print(f"Warning: {e}")
                except nx.NodeNotFound as e:
                    if verbose:
                        print(f"Warning: {e}")

        self.subG, self.connected_targets, self.shortest_paths_res = self.get_subnetwork(self.shortest_paths_res, sign_consistent=False)
        
        degrees = [deg for node, deg in self.subG.degree()]

        end_time = time.time()
        elapsed_time = end_time - start_time

        maxlength_connectedtargets = len(set(itertools.chain.from_iterable(self.connected_targets.values())))

        runinfo_entry = pd.Series({'analysis': 'shortest_paths', 'label': self.label, 'threshold': self.threshold, 'num_nodes': len(set(self.subG.nodes())), 'num_edges': len(set(self.subG.edges())), 'degrees': degrees, 'elapsed_time': elapsed_time, 'targets_connected': maxlength_connectedtargets}).to_frame().transpose()
        self.runinfo_df = pd.concat([self.runinfo_df, runinfo_entry], ignore_index=True)

        return self.shortest_paths_res

    def compute_paths(self, args):
        G, source, targets, cutoff = args
        paths_for_source = []
        connected_targets_for_source = []
        for target in targets:
            paths = list(nx.all_simple_paths(G, source=source, target=target, cutoff=cutoff))
            paths_for_source.extend(paths)
            if paths:
                connected_targets_for_source.append(target)
        return paths_for_source, connected_targets_for_source


    def all_paths(self, cutoff=None, verbose=False, num_processes=None):
        """
        Calculate all paths between sources and targets.

        Args:
            cutoff (int, optional): Cutoff for path length. If None, there's no cutoff.
            verbose (bool): If True, print warnings when no path is found to a given target.

        Returns:
            list: A list containing all paths.
        """
        self.label = f'{self.label}__allpaths'
        start_time = time.time()
        self.all_paths_res = []
        self.connected_all_path_targets = {}
        sources = list(self.source_dict.keys())
        targets = list(self.target_dict.keys())

        with Pool(processes=num_processes) as pool:
            all_args = [(self.subG, source, targets, cutoff) for source in sources]
            results = pool.map(self.compute_paths, all_args)

        for i, source in enumerate(sources):
            paths_for_source, connected_targets_for_source = results[i]
            self.all_paths_res.extend(paths_for_source)
            self.connected_all_path_targets[source] = connected_targets_for_source

        self.subG, self.connected_all_path_targets, self.all_paths_res = self.get_subnetwork(self.all_paths_res, sign_consistent=False)
        degrees = [deg for node, deg in self.subG.degree()]

        end_time = time.time()
        elapsed_time = end_time - start_time

        maxlength_connectedtargets = len(set(itertools.chain.from_iterable(self.connected_all_path_targets.values())))

        runinfo_entry = pd.Series({'analysis': 'all_paths', 'label': self.label, 'threshold': self.threshold, 'num_nodes': len(self.subG.nodes()), 'num_edges': len(self.subG.edges()), 'degrees': degrees, 'elapsed_time': elapsed_time, 'targets_connected': maxlength_connectedtargets}).to_frame().transpose()
        self.runinfo_df = pd.concat([self.runinfo_df, runinfo_entry], ignore_index=True)

        return self.all_paths_res

    

    def sign_consistency_check(self, paths):
        """
        Check the sign consistency of the shortest paths.

        Args:
            paths (list): A list of paths to check for sign consistency.

        Returns:
            list: A list containing the sign consistent paths.
        """
        self.label = f'{self.label}__sc'
        start_time = time.time()

        self.subG, self.connected_sc_targets, self.sc_paths_res = self.get_subnetwork(paths, sign_consistent=True)

        if not self.connected_sc_targets:
            maxlength_connectedtargets = 0
        else:
            maxlength_connectedtargets = len(set(itertools.chain.from_iterable(self.connected_sc_targets.values())))
        
        degrees = [deg for node, deg in self.subG.degree()]
        end_time = time.time()
        elapsed_time = end_time - start_time
        runinfo_entry = pd.Series({'analysis': 'sc_check', 'label': self.label, 'threshold': self.threshold, 'num_nodes': len(self.subG.nodes()), 'num_edges': len(self.subG.edges()), 'degrees': degrees, 'elapsed_time': elapsed_time, 'targets_connected': maxlength_connectedtargets}).to_frame().transpose()
        self.runinfo_df = pd.concat([self.runinfo_df, runinfo_entry], ignore_index=True)

        return self.sc_paths_res



    def reverse_graph(self, G):
        """
        Reverse the direction of all edges in the graph.
        """
        G = G.reverse()
        self.is_reversed = not self.is_reversed
        return G



    def to_SIFfile(self, paths, title):
        """
        Converts paths to SIF (Simple Interaction Format) and save as a file.

        Args:
            paths (list): List of paths to be converted.
            title (str): Title for the output file.

        Returns:
            DataFrame: A dataframe containing the SIF representation.
        """
        os.makedirs('./results', exist_ok=True)
        sif_tuples = []
        for path in paths:
            for i in range(len(path) - 1):
                interaction_type = 'P' if self.subG[path[i]][path[i+1]]['sign'] > 0 else 'N'
                sif_tuples.append((path[i], interaction_type, path[i+1]))

        sif_df = pd.DataFrame(sif_tuples, columns=['source', 'interaction', 'target']).drop_duplicates()
        sif_df.to_csv(title, sep='\t', index=None)
        return(sif_df)



    def visualize_graph(self, paths, title="Graph", is_sign_consistent=True):
        """
        Visualizes the graph via graphviz using the computed paths.

        Args:
            paths (list): List of paths to be visualized.
            title (str, optional): Filename for the visualization. Defaults to "Graph".
            is_sign_consistent (bool, optional): If True, only visualize sign consistent paths. Defaults to True.

        Returns:
            None
        """
        if len(paths) == 0:
            print('There were no sign consistent paths for the given perturbations and downstream effects.')
            return

        visualizer = GraphVisualizer(self)

        visualizer.visualize_graph(paths, title=title, is_sign_consistent=is_sign_consistent)



    def reachability_filter(self, G):
        """
        Filters the graph for reachability.

        Returns:
            None
        """
        reachable_nodes = set(self.source_dict.keys())
        for source in self.source_dict.keys():
            reachable_nodes.update(nx.descendants(G, source))
        
        subG = G.subgraph(reachable_nodes)

        self.label = f'{self.label}__reachability'

        return subG



    def network_batchrun(self, iter, cutoff=3, initial_threshold=0.01):
        """
        Executes a batch run for network analysis based on varying pagerank thresholds.
        Two execution threads are supported so far:
            1) Pagerank Thresholding > compute all paths > sign consistency check > shortest paths
            2) Pagerank Thresholding > shortest paths > sign consistency check

        Args:
            cutoff (int): The maximum depth to search paths in the network.
            initial_threshold (float): The starting pagerank threshold value.
            verbose (bool): Whether or not to print details during execution.

        Returns:
            None
        """

        
        self.label = f'{self.study_id}__{iter}'
        self.iter = iter
        self.subG = self.reachability_filter(self.G)
        
        try:
            self.pagerank_solver(personalize_for='source')
            self.pagerank_solver(personalize_for='target')
        except ZeroDivisionError as e:
            print(f"Warning: The network is too small. {e}")
            return
        initial_subG = self.subG
        initial_label = self.label
    
        self.threshold = initial_threshold

        while self.threshold >= 0:
            print('Computing path 1 with threshold', self.threshold)
            self.label = initial_label
            self.subG = initial_subG
            self.compute_overlap()
            shortest_paths = self.shortest_paths()
            self.to_SIFfile(shortest_paths, title=f'./results/{self.label}.sif')
            self.visualize_graph(shortest_paths, title=self.label, is_sign_consistent=True)
            shortest_sc_paths = self.sign_consistency_check(shortest_paths)
            self.to_SIFfile(shortest_sc_paths, title=f'./results/{self.label}.sif')
            self.visualize_graph(shortest_sc_paths, title=self.label, is_sign_consistent=True)
            self.threshold = round(self.threshold - 0.001, 3)
        
        self.threshold = initial_threshold

        while self.threshold > 0:
            print('Computing path 2 with threshold', self.threshold)
            self.label = initial_label
            self.subG = initial_subG
            self.compute_overlap()
            self.subG = self.reachability_filter(self.subG)
            all_paths = self.all_paths(cutoff=cutoff)
            self.to_SIFfile(all_paths, title=f'./results/{self.label}.sif')
            self.visualize_graph(all_paths, title=self.label, is_sign_consistent=True)
            all_sc_paths = self.sign_consistency_check(all_paths)
            self.to_SIFfile(all_sc_paths, title=f'./results/{self.label}.sif')
            self.visualize_graph(all_sc_paths, title=self.label, is_sign_consistent=True)
            self.threshold = round(self.threshold - 0.001, 3)

        self.runinfo_df.to_csv(f'./results/{self.study_id}__runinfo.csv', index=None)

        visualizer = GraphVisualizer(self)
        # visualizer.visualize_qc_thresholds()
        # visualizer.visualize_threshold_elbowplot()
        # visualizer.visualize_degrees()
        # visualizer.visualize_intersection()






class GraphVisualizer:
    """
    Incorporates different plots to compare the networks.
    """

    def __init__(self, graph_solver):
        """
        Initialize the graph visualizer with a graph solver.

        Args:
            graph_solver: An instance of a graph solver.
        """
        self.G = graph_solver.G
        self.subG = graph_solver.subG
        self.sources = list(graph_solver.source_dict.keys())
        self.targets = list(graph_solver.target_dict.keys())
        self.source_dict = graph_solver.source_dict
        self.target_dict = graph_solver.target_dict
        self.runinfo_df = graph_solver.runinfo_df
        self.study_id = graph_solver.study_id
        self.selected_thresholds = []
        self.iter = graph_solver.iter



    def visualize_graph(self, paths, title="Graph", is_sign_consistent=True):
        """
        Visualize the graph using the provided paths.

        Args:
            paths (list): Paths to be visualized.
            title (str, optional): Filename for the visualization. Defaults to "Graph".
            is_sign_consistent (bool, optional): If True, only visualize sign consistent paths. Defaults to True.
        """
        os.makedirs('./results', exist_ok=True)

        V = nx.DiGraph()

        for path in paths:
            for node in path:
                V.add_node(node)
                if len(path) > 1:
                    for i in range(len(path) - 1):
                        edge_data = self.subG.get_edge_data(path[i], path[i + 1])
                        V.add_edge(path[i], path[i + 1], **edge_data)

        sources = [s for s in self.sources if s in V.nodes()]
        targets = [path[-1] for path in paths]

        if len(V.edges()) > 600:
            print(f'The graph is too large to visualize. It has {len(V.edges())} edges.')
            return

        A = nx.nx_agraph.to_agraph(V)
        A.graph_attr['ratio'] = '0.70707'

        # Add an intermediary invisible node and edges for layout control
        with A.add_subgraph(name='cluster_sources', rank='min') as c:
            c.graph_attr['color'] = 'none'
            for s in sources:
                c.add_node(s)

        with A.add_subgraph(name='cluster_targets', rank='max') as c:
            c.graph_attr['color'] = 'none'
            for t in targets:
                c.add_node(t)

        for node in A.nodes():
            n = node.get_name()
            if n in sources:
                color = 'blue'
                if self.source_dict.get(n, 1) > 0 and is_sign_consistent:
                    fillcolor = 'lightblue'
                elif self.source_dict.get(n, 1) < 0 and is_sign_consistent:
                    fillcolor = 'coral'
                else:
                    fillcolor = 'white'
                node.attr['shape'] = 'box'
                node.attr['color'] = color
                node.attr['style'] = 'filled'
                node.attr['fillcolor'] = fillcolor
            elif n in targets:
                color = 'purple'
                if self.target_dict.get(n, 1) > 0 and is_sign_consistent:
                    fillcolor = 'lightblue'
                elif self.target_dict.get(n, 1) < 0 and is_sign_consistent:
                    fillcolor = 'coral'
                else:
                    fillcolor = 'white'
                node.attr['shape'] = 'ellipse'
                node.attr['color'] = color
                node.attr['style'] = 'filled'
                node.attr['fillcolor'] = fillcolor
            else:
                node.attr['shape'] = 'diamond'
                node.attr['color'] = 'gray'
                node.attr['style'] = 'filled'
                node.attr['fillcolor'] = 'white'
        
        for edge in A.edges():
            u, v = edge
            edge_data = V.get_edge_data(u, v)
            edge_color = 'green' if edge_data['sign'] == 1 else 'red'
            edge.attr['color'] = edge_color
        
        A.layout(prog='dot')
        file_path = f'./results/{title}.pdf'
        A.draw(file_path, prog='dot')



    def visualize_qc_thresholds(self):
        """
        Visualize number of nodes, edges, reached targets, and computational time over PageRank thresholds with grouped bars.
        """
        
        # Define labels for the bars
        network_labels = [
            ("__pagerank_{value}__shortest", "Shortest Path"), 
            ("__pagerank_{value}__shortest__sc", "Shortest Path + SC"), 
            ("__pagerank_{value}__reachability__allpaths", "All Paths"), 
            ("__pagerank_{value}__reachability__allpaths__sc", "All Paths + SC")
        ]
        
        thresholds = self.runinfo_df['threshold'].unique()
        
        fig, axs = plt.subplots(2, 2, figsize=(15, 10))
        
        bar_width = 0.15
        index = np.arange(len(thresholds))
        
        def plot_data(ax, column_name, title):
            for idx, (label, legend_label) in enumerate(network_labels):
                heights = []
                for threshold in thresholds:
                    current_label = label.replace("{value}", str(threshold))
                    subset = self.runinfo_df[self.runinfo_df['label'].str.contains(current_label)]
                    heights.append(subset[column_name].values[0] if not subset.empty else 0)
                
                bars = ax.bar(index + idx * bar_width, heights, bar_width, label=legend_label)
        
            ax.set_title(title)
            ax.set_xticks(index + bar_width * 1.5)
            ax.set_xticklabels([str(t) for t in thresholds])
            ax.legend()
        
        # Plot number of nodes
        plot_data(axs[0, 0], 'num_nodes', 'Number of Nodes')
        axs[0, 0].set_ylabel('Count')
        
        # Plot number of edges
        plot_data(axs[0, 1], 'num_edges', 'Number of Edges')
        axs[0, 1].set_ylabel('Count')
        
        # Plot number of reached targets
        plot_data(axs[1, 0], 'targets_connected', 'Number of Reached Targets')
        axs[1, 0].set_ylabel('Count')
        
        # Plot computational time
        plot_data(axs[1, 1], 'elapsed_time', 'Computational Time')
        axs[1, 1].set_ylabel('Time (s)')
        
        plt.tight_layout()
        plt.show()




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
        
        # Calculate unit vector of the line
        unit_line = (p2 - p1) / np.linalg.norm(p2 - p1)
        
        # Find the elbow by computing distance from each point to the line
        distances = []
        for i in range(len(x)):
            point = np.array([x[i], y[i]])
            distance = np.linalg.norm(np.cross(p2-p1, p1-point))/np.linalg.norm(p2-p1)
            distances.append(distance)
            
        # Return the index of the point with max distance which is the elbow
        return np.argmax(distances)



    def visualize_threshold_elbowplot(self):
        """
        Visualizes the threshold elbow plot for the graph with adjustments.
        The elbow in the plot indicates the optimal trade-off between the graph's 
        size and the number of connected targets.
        """
        
        network_labels = [
            ("__pagerank_{value}__shortest", "Shortest Path"), 
            ("__pagerank_{value}__shortest__sc", "Shortest Path + SC"), 
            ("__pagerank_{value}__reachability__allpaths", "All Paths"), 
            ("__pagerank_{value}__reachability__allpaths__sc", "All Paths + SC")
        ]
        
        thresholds_list = list(self.runinfo_df['threshold'].unique())
        selected_thresholds = [
            "__pagerank_0.0__shortest", 
            "__pagerank_0.0__shortest__sc"]  # To store threshold = 0
        
        fig, ax = plt.subplots(figsize=(10, 5))
        
        def plot_data(label, legend_label, color):
            num_edges = []
            missing_targets = []
            for threshold in thresholds_list:
                current_label = label.replace("{value}", str(threshold))
                subset = self.runinfo_df[self.runinfo_df['label'].str.contains(current_label)]
                
                if threshold == 0 and "allpaths" in current_label:
                    continue
                
                num_edges.append(subset['num_edges'].values[0] if not subset.empty else 0)
                missing_targets.append((1 - subset['targets_connected'].values[0] / len(self.target_dict)) * 100 if not subset.empty else 0)
            
            num_edges = [(e - min(num_edges)) / (max(num_edges) - min(num_edges)) for e in num_edges]
            
            elbow = self.find_elbow(num_edges, missing_targets)
            ax.plot(num_edges, missing_targets, '-o', label=legend_label, color=color)
            ax.scatter(num_edges[elbow], missing_targets[elbow], color='red', s=100)
            for idx, t in enumerate(thresholds_list):
                if t == 0 and "allpaths" in current_label:
                    continue
                ax.text(num_edges[idx], missing_targets[idx], str(t), fontsize=8)
            
            selected_thresholds.append(label.replace("{value}", str(thresholds_list[elbow])))

        colors = ['cornflowerblue', 'lightskyblue', 'lightgreen', 'plum']
        for idx, (label, legend_label) in enumerate(network_labels):
            plot_data(label, legend_label, colors[idx])

        ax.set_ylim(0, 100)
        ax.set_title('Rescaled Number of Edges vs Missing Targets, Elbow Plot')
        ax.set_ylabel('% Missing Targets')
        ax.set_xlabel('Rescaled Number of Edges (Min-Max Scaling)')
        ax.legend()
        
        self.selected_thresholds = selected_thresholds

        plt.show()



    def visualize_degrees(self):
        """
        Visualize the normalized degree distribution of the graph for selected thresholds.
        Provides statistical testing comparing the degree distribution of the 
        subnetworks.
        """
        
        formatted_data = []
        for threshold in self.selected_thresholds:
            subset = self.runinfo_df[self.runinfo_df['label'].str.contains(threshold)]
            if not subset.empty:
                degrees = subset['degrees'].iloc[0]
                total_nodes = len(degrees)  # Assuming 'degrees' contains a degree for each node in the subset
                for degree in degrees:
                    normalized_degree = degree / (2 * (total_nodes - 1))
                    formatted_data.append([threshold, normalized_degree])

        formatted_df = pd.DataFrame(formatted_data, columns=['Label', 'Normalized Degree'])

        plt.figure(figsize=(15, 10))
        ax = pt.RainCloud(data = formatted_df, x = 'Label', y = 'Normalized Degree', palette = "Set2",
                            orient = 'v', width_box= .1, width_viol=1)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45,  ha='right')
        test_results = st.add_stat_annotation(ax, data=formatted_df, x='Label', y='Normalized Degree', 
                                                box_pairs=[(self.selected_thresholds[i], self.selected_thresholds[j]) for i in range(len(self.selected_thresholds)) for j in range(i+1, len(self.selected_thresholds))],
                                                test='Mann-Whitney', text_format='star', loc='outside', verbose=2)

        plt.show()
    
    

    def get_threshold(self, filename):
        """
        Extracts the threshold value from a given filename.

        Args:
            filename (str): Name of the file.

        Returns:
            str: Extracted threshold value.
        """
        return filename.replace('.sif', '')



    def get_intersection_for_thresholds(self, thresholds, edges_dict):
        """
        Gets the intersection of edges for given thresholds.

        Args:
            thresholds (list): List of thresholds.
            edges_dict (dict): Dictionary containing edges for each threshold.

        Returns:
            int: Number of intersecting edges.
        """
        edge_sets = [edges_dict[thresh] for thresh in thresholds]
        return len(set.intersection(*edge_sets))



    def visualize_intersection(self):
        """
        Visualizes an upset plot showing the intersection of edges for graphs 
        built with the selected thresholds.
        """
        import os
        import itertools
        from upsetplot import plot
        
        directory = './results/'
        
        expected_files = [f'{self.study_id}__{self.iter}__reachability{threshold}.sif' 
                        for threshold in self.selected_thresholds]
        
        all_files = [f for f in os.listdir(directory) if f in expected_files]
        print(all_files)

        edges_dict = {}
        for file in all_files:
            path = os.path.join(directory, file)
            df = pd.read_csv(path, sep='\t')
            
            edges = set(tuple(row) for index, row in df.iterrows())
            edges_dict[self.get_threshold(file)] = edges

        intersection_dict = {}
        for r in range(1, len(edges_dict.keys()) + 1):
            for combination in itertools.combinations(edges_dict.keys(), r):
                index_values = [threshold in combination for threshold in edges_dict.keys()]
                intersection_dict[tuple(index_values)] = self.get_intersection_for_thresholds(combination, edges_dict)

        intersection_series = pd.Series(intersection_dict)
        intersection_series = intersection_series[intersection_series > 0]
        intersection_series.index.names = combination

        plot(intersection_series, sort_by='cardinality', show_counts=True, show_percentages=True)
        plt.show()