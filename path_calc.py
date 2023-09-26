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
            product_sign = 1
            source = path[0]
            target = path[-1]
            source_sign = self.source_dict[source]
            target_sign = self.target_dict[target]

            if source not in connected_targets:
                connected_targets[source] = []

            for i in range(len(path) - 1):
                edge_sign = self.subG.get_edge_data(path[i], path[i + 1])['sign']
                
                product_sign *= edge_sign

            if (sign_consistent and product_sign == source_sign * target_sign) or not sign_consistent:
                paths_res.append(path)
                connected_targets[source].append(target) if target not in connected_targets[source] else None
                for i in range(len(path) - 1):
                    edge_data = self.subG.get_edge_data(path[i], path[i + 1])
                    V.add_edge(path[i], path[i + 1], **edge_data)
        
        return V, connected_targets, paths_res



    def pagerank_solver(self, alpha=0.85, max_iter=100, tol=1.0e-6, nstart=None, weight='weight', dangling=None, personalize_for="source"):
        """
        Compute the PageRank values for nodes in the graph.

        Args:
            alpha (float): Damping factor for the PageRank algorithm.
            max_iter (int): Maximum number of iterations.
            tol (float): Tolerance to determine convergence.
            nstart (dict): Starting value of PageRank iteration for all nodes.
            weight (str): Edge data key to use as weight.
            dangling (Any): Value to use for dangling nodes.
            personalize_for (str): Personalize the PageRank by setting initial probabilities for either sources or targets.
        """
        if personalize_for == "source":
            personalized_prob = {n: 1/len(self.source_dict) for n in self.source_dict}
        elif personalize_for == "target":
            personalized_prob = {n: 1/len(self.target_dict) for n in self.target_dict}
        else:
            raise ValueError("personalize_for should be either 'source' or 'target'")
        
        if personalize_for == "target" and not self.is_reversed:
            self.reverse_graph()
        elif personalize_for == "source" and self.is_reversed:
            self.reverse_graph()
        
        personalized_prob.update({n: 0 for n in self.G.nodes() if n not in personalized_prob})
        
        pagerank = nx.pagerank(self.G, alpha=alpha, max_iter=max_iter, personalization=personalized_prob, tol=tol, nstart=nstart, weight=weight, dangling=dangling)
        
        self.pagerank_values[personalize_for] = pagerank
        
        for node, pr_value in pagerank.items():
            attribute_name = 'pagerank_from_targets' if personalize_for == "target" else 'pagerank_from_sources'
            self.G.nodes[node][attribute_name] = pr_value
        
        if personalize_for == "target" and self.is_reversed:
            self.reverse_graph()



    def compute_overlap(self):
        """
        Compute the overlap of nodes that exceed the PageRank threshold from sources and targets.

        Returns:
            tuple: Contains nodes above threshold from sources, nodes above threshold from targets, and overlapping nodes.
        """
        nodes_above_threshold_from_sources = {node for node, data in self.G.nodes(data=True) if data.get('pagerank_from_sources') > self.threshold}
        nodes_above_threshold_from_targets = {node for node, data in self.G.nodes(data=True) if data.get('pagerank_from_targets') > self.threshold}
        overlap = nodes_above_threshold_from_sources.intersection(nodes_above_threshold_from_targets)
        nodes_to_include = nodes_above_threshold_from_sources.union(nodes_above_threshold_from_targets)
        self.subG = self.G.subgraph(nodes_to_include)
        self.label = f'{self.label}__{self.threshold}'
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



    def all_paths(self, cutoff=None, verbose=False):
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
        for source_node in sources:
            if source_node not in self.connected_all_path_targets:
                self.connected_all_path_targets[source_node] = []

            for target_node in targets:
                try:
                    paths = list(nx.all_simple_paths(self.subG, source=source_node, target=target_node, cutoff=cutoff))
                    self.all_paths_res.extend(paths)
                    if paths:
                        self.connected_all_path_targets[source_node].append(target_node)
                except Exception as e:
                    if verbose:
                        print(f"Warning: {e}")

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



    def reverse_graph(self):
        """
        Reverse the direction of all edges in the graph.
        """
        self.G = self.G.reverse()
        self.is_reversed = not self.is_reversed 



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



    def visualize_qcplots(self):
        """
        Visualizes quality control plots for the graph.

        Returns:
            None
        """
        visualizer = GraphVisualizer(self)

        visualizer.visualize_size_thresholds()
        visualizer.visualize_threshold_elbowplot()
        visualizer.visualize_comptime()



    def visualize_degrees(self):
        """
        Visualizes the degree distribution of the graph for selected thresholds.

        Returns:
            None
        """
        visualizer = GraphVisualizer(self)
        visualizer.visualize_degrees()



    def visualize_intersection(self):
        """
        Visualizes the intersection of edges for selected thresholds.

        Returns:
            None
        """
        visualizer = GraphVisualizer(self)
        visualizer.visualize_intersection()



    def network_batchrun(self, cutoff=3, initial_threshold=0.01, verbose=False):
        """
        Executes a batch run for network analysis based on varying thresholds.

        Args:
            cutoff (int): The maximum depth to search paths in the network.
            initial_threshold (float): The starting pagerank threshold value.
            verbose (bool): Whether or not to print details during execution.

        Returns:
            None
        """
        self.pagerank_solver(personalize_for='source')
        self.pagerank_solver(personalize_for='target')
    
        self.threshold = initial_threshold
        while self.threshold >=0:
            self.label = 'pagerank'
            self.compute_overlap()
            self.all_paths(cutoff=cutoff)
            self.sign_consistency_check(self.all_paths_res)
            paths = self.shortest_paths(verbose=verbose)
            self.to_SIFfile(paths, title=f'./results/{self.study_id}__{self.label}.sif')
            self.visualize_graph(paths, title=self.label, is_sign_consistent=True)
            self.threshold = round(self.threshold - 0.001, 3)
        
        self.threshold = initial_threshold
        while self.threshold >=0:
            self.label = 'pagerank'
            self.compute_overlap()
            self.shortest_paths(verbose=verbose)
            paths = self.sign_consistency_check(self.shortest_paths_res)
            self.to_SIFfile(paths, title=f'./results/{self.study_id}__{self.label}.sif')
            self.visualize_graph(paths, title=self.label, is_sign_consistent=True)
            self.threshold = round(self.threshold - 0.001, 3)

        visualizer = GraphVisualizer(self)
        visualizer.visualize_size_thresholds()
        visualizer.visualize_threshold_elbowplot()
        visualizer.visualize_comptime()
        visualizer.visualize_degrees()
        visualizer.visualize_intersection()





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
        file_path = f'./results/{self.study_id}_{title}.pdf'
        A.draw(file_path, prog='dot')



    def visualize_size_thresholds(self):
        """
        Visualize number of nodes, edges, and % connected targets over PageRank thresholds.
        """
        thresholds_allpaths = self.runinfo_df['threshold'][self.runinfo_df['label'].str.endswith('allpaths__sc__shortest')].tolist()
        num_nodes_allpaths = self.runinfo_df['num_nodes'][self.runinfo_df['label'].str.endswith('allpaths__sc__shortest')].tolist()
        num_edges_allpaths = self.runinfo_df['num_edges'][self.runinfo_df['label'].str.endswith('allpaths__sc__shortest')].tolist()
        num_targets_allpaths = self.runinfo_df['targets_connected'][self.runinfo_df['label'].str.endswith('allpaths__sc__shortest')].tolist()

        thresholds_shortestpaths = self.runinfo_df['threshold'][self.runinfo_df['label'].str.endswith('shortest__sc')].tolist()
        num_nodes_shortestpaths = self.runinfo_df['num_nodes'][self.runinfo_df['label'].str.endswith('shortest__sc')].tolist()
        num_edges_shortestpaths = self.runinfo_df['num_edges'][self.runinfo_df['label'].str.endswith('shortest__sc')].tolist()
        num_targets_shortestpaths = self.runinfo_df['targets_connected'][self.runinfo_df['label'].str.endswith('shortest__sc')].tolist()

        perc_connected_targets_allpaths = [round((num_targets_allpaths[i]/len(self.target_dict))*100, 2) for i in range(len(num_targets_allpaths))]
        perc_connected_targets_shortestpaths = [round((num_targets_shortestpaths[i]/len(self.target_dict))*100, 2) for i in range(len(num_targets_shortestpaths))]

        plt.figure(figsize=(10, 5))

        # Primary y-axis (on the left)
        ax1 = plt.gca()  # get the current axis
        ax1.plot(thresholds_allpaths, num_nodes_allpaths, '-o', label="Number of Nodes, all paths")
        ax1.plot(thresholds_allpaths, num_edges_allpaths, '-o', label="Number of Edges, all paths")
        ax1.plot(thresholds_shortestpaths, num_nodes_shortestpaths, '-o', label="Number of Nodes, shortest paths")
        ax1.plot(thresholds_shortestpaths, num_edges_shortestpaths, '-o', label="Number of Edges, shortest paths")
        ax1.set_xlabel('Threshold')
        ax1.set_ylabel('Count')
        ax1.tick_params(axis='y')
        ax1.legend(loc='right')
        ax1.set_ylim(0)

        # Secondary y-axis (on the right)
        ax2 = ax1.twinx()  # Create a twin y-axis sharing the same x-axis
        ax2.plot(thresholds_allpaths, perc_connected_targets_allpaths, '-o', label="Number of Targets connected, all paths", color = 'olive')
        ax2.plot(thresholds_shortestpaths, perc_connected_targets_shortestpaths, '-o', label="Number of Targets connected, shortest paths", color = 'purple')
        ax2.set_ylabel('% connected Targets')
        ax2.set_ylim(0, 100)
        ax2.tick_params(axis='y')
        ax2.legend(loc='upper right')

        plt.title('PageRank: Number of Nodes, Edges, and % connected Targets over Thresholds')
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
        Visualizes the threshold elbow plot for the graph.

        The elbow in the plot indicates the optimal trade-off between the graph's 
        size and the number of connected targets.
        """
        thresholds = self.runinfo_df['threshold'][self.runinfo_df['label'].str.endswith('allpaths__sc__shortest')].tolist()
        num_edges_allpaths = self.runinfo_df['num_edges'][self.runinfo_df['label'].str.endswith('allpaths__sc__shortest')].tolist()
        num_targets_allpaths = self.runinfo_df['targets_connected'][self.runinfo_df['label'].str.endswith('allpaths__sc__shortest')].tolist()
        num_edges_shortestpaths = self.runinfo_df['num_edges'][self.runinfo_df['label'].str.endswith('shortest__sc')].tolist()
        num_targets_shortestpaths = self.runinfo_df['targets_connected'][self.runinfo_df['label'].str.endswith('shortest__sc')].tolist()

        # Find elbows
        elbow_allpaths = self.find_elbow(num_edges_allpaths, [(1-t/len(self.target_dict))*100 for t in num_targets_allpaths])
        elbow_shortestpaths = self.find_elbow(num_edges_shortestpaths, [(1-t/len(self.target_dict))*100 for t in num_targets_shortestpaths])

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
        missing_targets_allpaths = [(1-t/len(self.target_dict))*100 for t in num_targets_allpaths]
        missing_targets_shortestpaths = [(1-t/len(self.target_dict))*100 for t in num_targets_shortestpaths]
        ax1.plot(num_edges_allpaths, missing_targets_allpaths, '-o')
        ax2.plot(num_edges_shortestpaths, missing_targets_shortestpaths, '-o')
        
        # Highlight elbows
        ax1.scatter(num_edges_allpaths[elbow_allpaths], missing_targets_allpaths[elbow_allpaths], color='red', s=100)
        ax2.scatter(num_edges_shortestpaths[elbow_shortestpaths], missing_targets_shortestpaths[elbow_shortestpaths], color='red', s=100)

        for t in thresholds:
            ax1.text(num_edges_allpaths[thresholds.index(t)], missing_targets_allpaths[thresholds.index(t)], str(t))
            ax2.text(num_edges_shortestpaths[thresholds.index(t)], missing_targets_shortestpaths[thresholds.index(t)], str(t))
        ax1.set_ylim(0, 100)
        ax2.set_ylim(0, 100)

        ax1.set_title('Number of Edges vs missing Targets, Elbow Plot, All paths')
        ax2.set_title('Number of Edges vs missing Targets, Elbow Plot, Shortest paths')
        ax1.set_ylabel('% missing Targets')
        ax1.set_xlabel('Number of Edges')
        ax2.set_ylabel('% missing Targets')
        ax2.set_xlabel('Number of Edges')

        label_allpaths = self.runinfo_df['label'][(self.runinfo_df['threshold'] == thresholds[elbow_allpaths]) & 
                                            (self.runinfo_df['label'].str.endswith('allpaths__sc__shortest'))].iloc[0]
        label_shortestpaths = self.runinfo_df['label'][(self.runinfo_df['threshold'] ==  thresholds[elbow_allpaths]) & 
                                                (self.runinfo_df['label'].str.endswith('shortest__sc'))].iloc[0]
        label_p0_allpaths = self.runinfo_df['label'][(self.runinfo_df['threshold'] == 0) & 
                                            (self.runinfo_df['label'].str.endswith('allpaths__sc__shortest'))].iloc[0]
        label_p0_shortestpaths = self.runinfo_df['label'][(self.runinfo_df['threshold'] ==  0) & 
                                                (self.runinfo_df['label'].str.endswith('shortest__sc'))].iloc[0]

        self.selected_thresholds = [label_allpaths, label_shortestpaths, label_p0_allpaths, label_p0_shortestpaths]

        plt.show()
    


    def visualize_comptime(self):
        """
        Visualize the computation time of the graph.
        """
        thresholds_allpaths = self.runinfo_df['threshold'][self.runinfo_df['label'].str.endswith('allpaths__sc__shortest')].tolist()
        thresholds_shortestpaths = self.runinfo_df['threshold'][self.runinfo_df['label'].str.endswith('shortest__sc')].tolist()

        computation_times_shortestpaths = [x + y for x, y in zip(self.runinfo_df['elapsed_time'][(self.runinfo_df['label'].str.endswith('shortest')) & (~self.runinfo_df['label'].str.contains('allpaths'))].tolist(), self.runinfo_df['elapsed_time'][self.runinfo_df['label'].str.endswith('shortest__sc')].tolist())]
        computation_times_allpaths = [x + y + z for x, y, z in zip(self.runinfo_df['elapsed_time'][self.runinfo_df['label'].str.endswith('allpaths__sc__shortest')].tolist(), self.runinfo_df['elapsed_time'][self.runinfo_df['label'].str.endswith('allpaths__sc')].tolist(), self.runinfo_df['elapsed_time'][self.runinfo_df['label'].str.endswith('allpaths')].tolist())]
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
        ax1.plot(thresholds_allpaths, computation_times_allpaths, '-o', color="red")
        ax2.plot(thresholds_shortestpaths, computation_times_shortestpaths, '-o', color="red")
        ax1.set_title('Computation time, All paths')
        ax2.set_title('Computation time, Shortest paths')
        ax1.set_ylabel('Time (s)')
        ax1.set_xlabel('Threshold')
        ax2.set_ylabel('Time (s)')
        ax2.set_xlabel('Threshold')
        fig.show()
        


    def visualize_degrees(self):
        """
        Visualize the degree distribution of the graph for selected thresholds.
        Provides statistical testing comparing the degree distribution of the 
        subnetworks.
        """
        formatted_data = []
        for threshold in self.selected_thresholds:
            subset = self.runinfo_df[self.runinfo_df['label'] == threshold]
            if not subset.empty:
                degrees = subset['degrees'].iloc[0]
                for degree in degrees:
                    formatted_data.append([threshold, degree])

        formatted_df = pd.DataFrame(formatted_data, columns=['Threshold', 'Degree'])

        plt.figure(figsize=(10, 5))
        ax = pt.RainCloud(data = formatted_df, x = 'Threshold', y = 'Degree', palette = "Set2",
                        orient = 'v', width_box= .1, width_viol=1)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        test_results = st.add_stat_annotation(ax, data=formatted_df, x='Threshold', y='Degree', 
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

        Args:
            selected_thresholds (list): List of thresholds to visualize.
        """
        directory = './results/'

        expected_files = [f'{self.study_id}__{threshold}.sif' for threshold in self.selected_thresholds]
        all_files = [f for f in os.listdir(directory) if f in expected_files]

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


        intersection_series = pd.Series(intersection_dict, )
        intersection_series = intersection_series[intersection_series > 0]
        intersection_series.index.names = combination

        plot(intersection_series, sort_by='cardinality', show_counts=True, show_percentages=True)
        plt.show()