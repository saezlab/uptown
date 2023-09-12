# shortest_path.py

import pandas as pd
import networkx as nx
import pygraphviz as pgv
import numpy as np
from matplotlib import pyplot as plt
import os
from upsetplot import plot
import itertools


class GraphSolverCore:
    """
    Base class for solving problems related to graphs.
    """
    def __init__(self, G, study_id = None):
        """
        Initialize the base graph solver with a given graph.

        :param G: A networkx graph object.
        """
        self.G = G
        self.study_id = study_id




class BaselineSolver(GraphSolverCore):
    """
    Solver for finding signed paths in a graph.
    """

    def __init__(self, G, study_id = None):
        """
        Initialize the signed path solver with a given graph.

        :param G: A networkx graph object.
        """

        super().__init__(G, study_id)
        self.source_dict = {}
        self.target_dict = {}
        self.shortest_paths_res = []
        self.shortest_sc_paths_res = []
        self.connected_targets = {}
        self.connected_sc_targets = {}
        self.runinfo_dict = {}



    def shortest_paths(self, label, verbose = False):
        """
        Calculate the shortest paths between sources and targets.

        :param label: A label for the current run.
        :param verbose: If True, print warnings when no path is found to a given target.
        :return: A tuple containing the shortest paths and run information.
        """
        self.shortest_paths_res = []
        self.connected_targets = {}
        sources = list(self.source_dict.keys())
        targets = list(self.target_dict.keys())
        for source_node in sources:
            if source_node not in self.connected_targets:
                self.connected_targets[source_node] = []

            for target_node in targets:
                try:
                    self.shortest_paths_res.extend([p for p in nx.all_shortest_paths(self.G, source=source_node, target=target_node, weight='weight')])
                    self.connected_targets[source_node].append(target_node)
                except nx.NetworkXNoPath as e:
                    if verbose:
                        print(f"Warning: {e}")

        
        degrees = [deg for node, deg in self.G.degree()]

        runinfo = {
            label: {
                'num_nodes': len(self.G.nodes()),
                'num_edges': len(self.G.edges()),
                'degrees': degrees
            }
        }

        return self.shortest_paths_res, runinfo
        


    def sign_consistency_check(self, label):
        """
        Check the sign consistency of the shortest paths.

        :param label: A label for the current run.
        :return: A tuple containing the sign consistent paths and run information.
        """
        self.shortest_sc_paths_res = []
        self.connected_sc_targets = {}
        V = nx.DiGraph()
        for path in self.shortest_paths_res:
            product_sign = 1
            source = path[0]
            target = path[-1]
            source_sign = self.source_dict[source]
            target_sign = self.target_dict[target]
        
            if source not in self.connected_sc_targets:
                self.connected_sc_targets[source] = []

            for i in range(len(path) - 1):
                edge_sign = self.G.get_edge_data(path[i], path[i + 1])['sign']
                edge_data = self.G.get_edge_data(path[i], path[i + 1])
                product_sign *= edge_sign
                V.add_edge(path[i], path[i + 1], **edge_data)
            
            if product_sign == source_sign * target_sign:
                self.shortest_sc_paths_res.append(path)
                self.connected_sc_targets[source].append(target)
            self.connected_sc_targets[source] = list(set(self.connected_sc_targets[source]))
        
        degrees = [deg for node, deg in V.degree()]
        
        runinfo = {
            label: {
                'num_nodes': len(V.nodes()),
                'num_edges': len(V.edges()),
                'degrees': degrees
            }
        }

        return self.shortest_sc_paths_res, runinfo
    

    
    def to_SIFfile(self, paths, title):
        """
        Convert paths to SIF (Simple Interaction Format) and save as a file.

        :param paths: List of paths to be converted.
        :param title: Title for the output file.
        :return: A dataframe containing the SIF representation.
        """
        sif_tuples = []
        for path in paths:
            for i in range(len(path) - 1):
                interaction_type = 'P' if self.G[path[i]][path[i+1]]['sign'] > 0 else 'N'
                sif_tuples.append((path[i], interaction_type, path[i+1]))

        sif_df = pd.DataFrame(sif_tuples, columns=['source', 'interaction', 'target']).drop_duplicates()
        sif_df.to_csv(title, sep='\t', index=None)
        return(sif_df)
    


    def visualize_graph(self, title="Graph", export_sif=False, is_sign_consistent=True):
        """
        Visualize the graph via graphviz using the computed paths.

        :param title: Filename for the visualization.
        :param export_sif: If True, export the visualization in SIF format.
        :param is_sign_consistent: If True, only visualize sign consistent paths.
        """
        if is_sign_consistent and len(self.shortest_sc_paths_res) == 0:
            print('There were no sign consistent paths for the given perturbations and downstream effects.')
            return
        
        paths = self.shortest_sc_paths_res if is_sign_consistent else self.shortest_paths_res

        visualizer = GraphVisualizer(self)

        visualizer.visualize_graph(paths, title=title, is_sign_consistent=is_sign_consistent)  

        self.to_SIFfile(paths, title=f'./results/{self.study_id}_{title}.sif') if export_sif else None




class GraphVisualizer:
    """
    Class that incorporates different plots to compare the networks.
    """

    def __init__(self, graph_solver):
        """
        Initialize the graph visualizer with a given graph solver.

        :param graph_solver: An instance of a graph solver.
        """
        self.G = graph_solver.G
        self.sources = list(graph_solver.source_dict.keys())
        self.targets = list(graph_solver.target_dict.keys())
        self.source_dict = graph_solver.source_dict
        self.target_dict = graph_solver.target_dict
        self.runinfo_dict = graph_solver.runinfo_dict
        self.study_id = graph_solver.study_id



    def visualize_graph(self, paths, title="Graph", is_sign_consistent=True):
        """
        Visualize the graph using the given paths.

        :param paths: List of paths to be visualized.
        :param title: Filename for the visualization.
        :param is_sign_consistent: If True, only visualize sign consistent paths.
        """
        os.makedirs('./results', exist_ok=True)

        V = nx.DiGraph()

        for path in paths:
            for node in path:
                V.add_node(node)
                if len(path) > 1:
                    for i in range(len(path) - 1):
                        edge_data = self.G.get_edge_data(path[i], path[i + 1])
                        V.add_edge(path[i], path[i + 1], **edge_data)

        sources = [s for s in self.sources if s in V.nodes()]
        targets = [path[-1] for path in paths]

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



    def visualize_size_thresholds(self, title = "SizeThresholds", is_sign_consistent=True):
        """
        Visualize the number of nodes, edges, and % connected targets over PageRank thresholds.

        :param title: Title for the visualization.
        :param is_sign_consistent: If True, only visualize sign consistent paths.
        """
        if not is_sign_consistent:
            thresholds = list({k: v for k, v in self.runinfo_dict.items() if k.startswith('pagerank_') and not k.startswith('pagerank_sc')})
            thresholds_float = [float(t.split('_')[1]) for t in thresholds]
        elif is_sign_consistent:
            thresholds = list({k: v for k, v in self.runinfo_dict.items() if k.startswith('pagerank_sc')})
            thresholds_float = [float(t.split('_')[2]) for t in thresholds]
        
        num_nodes = [self.runinfo_dict[t]['num_nodes'] for t in thresholds]
        num_edges = [self.runinfo_dict[t]['num_edges'] for t in thresholds]
        num_targets = [self.runinfo_dict[t]['targets_connected']/len(self.target_dict)*100 for t in thresholds]

        

        plt.figure(figsize=(10, 5))

        # Primary y-axis (on the left)
        ax1 = plt.gca()  # get the current axis
        ax1.plot(thresholds_float, num_nodes, '-o', label="Number of Nodes")
        ax1.plot(thresholds_float, num_edges, '-o', label="Number of Edges")
        ax1.set_xlabel('Threshold')
        ax1.set_ylabel('Count')
        ax1.tick_params(axis='y')
        ax1.legend(loc='upper center')

        # Secondary y-axis (on the right)
        ax2 = ax1.twinx()  # Create a twin y-axis sharing the same x-axis
        ax2.plot(thresholds_float, num_targets, '-o', label="Number of Targets connected", color='olive')
        ax2.set_ylabel('% connected Targets')
        ax2.set_ylim(0, 100)
        ax2.tick_params(axis='y')
        ax2.legend(loc='upper right')

        plt.title('PageRank: Number of Nodes, Edges, and % connected Targets over Thresholds')
        plt.show()



    def visualize_threshold_elbowplot(self, title = "Elbowplot", is_sign_consistent=True):
        """
        Visualize the threshold elbow plot of the graph. In other words, 
        the elbow indicates the best ratio size/number of connected targets.

        :param title: Title for the visualization.
        :param is_sign_consistent: If True, only visualize sign consistent paths.
        """
        if not is_sign_consistent:
            thresholds = list({k: v for k, v in self.runinfo_dict.items() if k.startswith('pagerank_') and not k.startswith('pagerank_sc')})
        elif is_sign_consistent:
            thresholds = list({k: v for k, v in self.runinfo_dict.items() if k.startswith('pagerank_sc')})

        num_edges = [self.runinfo_dict[t]['num_edges'] for t in thresholds]
        num_targets = [self.runinfo_dict[t]['targets_connected']/46*100 for t in thresholds]

        plt.figure(figsize=(10, 5))
        ax = plt.gca()
        missing_targets = [100-t for t in num_targets]
        ax.plot(num_edges, missing_targets, '-o')
        for t in thresholds:
            ax.text(num_edges[thresholds.index(t)], missing_targets[thresholds.index(t)], str(t))

        plt.title('Number of Edges vs missing Targets, Elbow Plot')
        plt.ylabel('% missing Targets')
        plt.xlabel('Number of Edges')
        plt.show()
    


    def visualize_comptime(self, is_sign_consistent=True):
        """
        Visualize the computation time of the graph.

        :param is_sign_consistent: If True, only visualize sign consistent paths.
        """
        if not is_sign_consistent:
            thresholds = list({k: v for k, v in self.runinfo_dict.items() if k.startswith('pagerank_') and not k.startswith('pagerank_sc')})
            thresholds_float = [float(t.split('_')[1]) for t in thresholds]
        elif is_sign_consistent:
            thresholds = list({k: v for k, v in self.runinfo_dict.items() if k.startswith('pagerank_sc')})
            thresholds_float = [float(t.split('_')[2]) for t in thresholds]
 
        computation_times = [self.runinfo_dict[t]['computation_time'] for t in thresholds]
        plt.figure(figsize=(10, 5))
        plt.plot(thresholds_float, computation_times, '-o', color="red")
        plt.xlabel('Threshold')
        plt.ylabel('Computation Time (s)')
        plt.title('Computation Time over Thresholds')
        plt.show()
        


    def visualize_degrees(self, selected_thresholds):
        """
        Visualize the degree distribution of the graph for selected thresholds.

        :param selected_thresholds: List of thresholds to visualize.
        """
        plt.figure(figsize=(10, 5))
        for threshold in selected_thresholds:
            data = self.runinfo_dict.get(threshold)
            if data:
                norm_degrees = np.ones_like(data['degrees']) / len(data['degrees'])
                plt.hist(data['degrees'], label=threshold, alpha=0.3, bins=np.linspace(0, 100, 101), weights=norm_degrees)

        plt.xlabel('Degree quantile')
        plt.ylabel('Frequency')
        plt.title('Degree Distribution across all Thresholds')
        plt.legend(loc='upper right')
        plt.show()
    
    

    def get_threshold(self, filename):
        """
        Extract the threshold value from a given filename.

        :param filename: Name of the file.
        :return: Extracted threshold value.
        """
        return filename.replace('.sif', '')



    def get_intersection_for_thresholds(self, thresholds, edges_dict):
        """
        Get the intersection of edges for given thresholds.

        :param thresholds: List of thresholds.
        :param edges_dict: Dictionary containing edges for each threshold.
        :return: Number of intersecting edges.
        """
        edge_sets = [edges_dict[thresh] for thresh in thresholds]
        return len(set.intersection(*edge_sets))



    def visualize_intersection(self, selected_thresholds):
        """
        Visualize an upset plot showing the intersection of edges 
        for graphs built with the selected thresholds.

        :param selected_thresholds: List of thresholds to visualize.
        """
        directory = './results/'

        expected_files = [f'{self.study_id}_{threshold}.sif' for threshold in selected_thresholds]
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