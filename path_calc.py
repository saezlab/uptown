import pandas as pd
import networkx as nx
import pygraphviz as pgv
import numpy as np
from matplotlib import pyplot as plt
import os
from upsetplot import plot
import itertools
import time
import ptitprince as pt


class Solver:
    """
    Solver for finding signed paths in a graph.
    """

    def __init__(self, G, study_id = None):
        """
        Initialize the signed path solver with a given graph.

        :param G: A networkx graph object.
        """

        self.G = G
        self.study_id = study_id
        self.label = 'pagerank'
        self.source_dict = {}
        self.target_dict = {}
        self.threshold = 0
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


    def pagerank_solver(self, alpha=0.85, max_iter=100, tol=1.0e-6, nstart=None, weight='weight', dangling=None, personalize_for="source"):
        """
        Compute the PageRank values for nodes in the graph.

        :param alpha: Damping factor for the PageRank algorithm.
        :param max_iter: Maximum number of iterations.
        :param tol: Tolerance to determine convergence.
        :param nstart: Starting value of PageRank iteration for all nodes.
        :param weight: Edge data key to use as weight.
        :param dangling: Value to use for dangling nodes.
        :param personalize_for: Personalize the PageRank by setting initial probabilities for either sources or targets.
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

        :return: A tuple containing nodes above threshold from sources, nodes above threshold from targets, and overlapping nodes.
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

        :param label: A label for the current run.
        :param verbose: If True, print warnings when no path is found to a given target.
        :return: A tuple containing the shortest paths and run information.
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

        
        degrees = [deg for node, deg in self.subG.degree()]

        end_time = time.time()
        elapsed_time = end_time - start_time

        maxlength_connectedtargets = max(len(targets) for targets in self.connected_targets.values())

        runinfo_entry = pd.Series({'analysis': 'shortest_paths', 'label': self.label, 'threshold': self.threshold, 'num_nodes': len(self.subG.nodes()), 'num_edges': len(self.subG.edges()), 'degrees': degrees, 'elapsed_time': elapsed_time, 'targets_connected': maxlength_connectedtargets})
        self.runinfo_df = self.runinfo_df.append(runinfo_entry, ignore_index=True)

        self.to_SIFfile(self.shortest_paths_res, title=f'./results/{self.study_id}__{self.label}.sif')

        return self.shortest_paths_res


    def all_paths(self, cutoff=None, verbose=False):
        """
        Calculate all paths between sources and targets.

        :param label: A label for the current run.
        :param cutoff: Cutoff for path length. If None, there's no cutoff.
        :param verbose: If True, print warnings when no path is found to a given target.
        :return: A tuple containing all paths and run information.
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

        degrees = [deg for node, deg in self.subG.degree()]

        end_time = time.time()
        elapsed_time = end_time - start_time

        maxlength_connectedtargets = max(len(targets) for targets in self.connected_all_path_targets.values())

        runinfo_entry = pd.Series({'analysis': 'all_paths', 'label': self.label, 'threshold': self.threshold, 'num_nodes': len(self.subG.nodes()), 'num_edges': len(self.subG.edges()), 'degrees': degrees, 'elapsed_time': elapsed_time, 'targets_connected': maxlength_connectedtargets})
        self.runinfo_df = self.runinfo_df.append(runinfo_entry, ignore_index=True)

        self.to_SIFfile(self.all_paths_res, title=f'./results/{self.study_id}__{self.label}.sif')

        return self.all_paths_res

    
    def sign_consistency_check(self, paths):
        """
        Check the sign consistency of the shortest paths.

        :param label: A label for the current run.
        :return: A tuple containing the sign consistent paths and run information.
        """
        self.label = f'{self.label}__sc'
        start_time = time.time()
        self.sc_paths_res = []
        self.connected_sc_targets = {}
        V = nx.DiGraph()
        for path in paths:
            product_sign = 1
            source = path[0]
            target = path[-1]
            source_sign = self.source_dict[source]
            target_sign = self.target_dict[target]

            if source not in self.connected_sc_targets:
                self.connected_sc_targets[source] = []

            for i in range(len(path) - 1):
                edge_sign = self.subG.get_edge_data(path[i], path[i + 1])['sign']
                edge_data = self.subG.get_edge_data(path[i], path[i + 1])
                product_sign *= edge_sign
                V.add_edge(path[i], path[i + 1], **edge_data)

            if product_sign == source_sign * target_sign:
                self.sc_paths_res.append(path)
                self.connected_sc_targets[source].append(target)
            self.connected_sc_targets[source] = list(set(self.connected_sc_targets[source]))

        self.subG = V

        maxlength_connectedtargets = max(len(targets) for targets in self.connected_sc_targets.values())
        
        degrees = [deg for node, deg in self.subG.degree()]
        end_time = time.time()
        elapsed_time = end_time - start_time
        runinfo_entry = pd.Series({'analysis': 'sc_check', 'label': self.label, 'threshold': self.threshold, 'num_nodes': len(self.subG.nodes()), 'num_edges': len(self.subG.edges()), 'degrees': degrees, 'elapsed_time': elapsed_time, 'targets_connected': maxlength_connectedtargets})
        self.runinfo_df = self.runinfo_df.append(runinfo_entry, ignore_index=True)
        
        self.to_SIFfile(self.sc_paths_res, title=f'./results/{self.study_id}__{self.label}.sif')

        return self.sc_paths_res


    def reverse_graph(self):
        """
        Reverse the direction of all edges in the graph.
        """
        self.G = self.G.reverse()
        self.is_reversed = not self.is_reversed 


    def to_SIFfile(self, paths, title):
        """
        Convert paths to SIF (Simple Interaction Format) and save as a file.

        :param paths: List of paths to be converted.
        :param title: Title for the output file.
        :return: A dataframe containing the SIF representation.
        """
        directory = './results/'
        if not os.path.exists(directory):
            os.makedirs(directory)
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
        Visualize the graph via graphviz using the computed paths.

        :param title: Filename for the visualization.
        :param export_sif: If True, export the visualization in SIF format.
        :param is_sign_consistent: If True, only visualize sign consistent paths.
        """
        if is_sign_consistent and len(self.sc_paths_res) == 0:
            print('There were no sign consistent paths for the given perturbations and downstream effects.')
            return

        visualizer = GraphVisualizer(self)

        visualizer.visualize_graph(paths, title=title, is_sign_consistent=is_sign_consistent)

    def visualize_qcplots(self):
        """
        Visualize quality control plots for the graph. Wrapper around the 
        visualize_size_thresholds, visualize_threshold_elbowplot and visualize_comptime 
        methods from GraphVisualizer class.

        :param title: Title for the visualization.
        :param is_sign_consistent: If True, only visualize sign consistent paths.
        """
        visualizer = GraphVisualizer(self)

        visualizer.visualize_size_thresholds()
        visualizer.visualize_threshold_elbowplot()
        visualizer.visualize_comptime()



    def visualize_degrees(self, selected_thresholds):
        """
        Visualize the degree distribution of the graph for selected thresholds. Wrapper 
        around the visualize_degrees method from GraphVisualizer class.

        :param selected_thresholds: List of thresholds to visualize.
        """
        visualizer = GraphVisualizer(self)
        visualizer.visualize_degrees(selected_thresholds)



    def visualize_intersection(self, selected_thresholds):
        """
        Visualize the intersection of edges for selected thresholds. Wrapper around the 
        visualize_intersection method from GraphVisualizer class.

        :param selected_thresholds: List of thresholds to visualize.
        """
        visualizer = GraphVisualizer(self)
        visualizer.visualize_intersection(selected_thresholds)




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
        self.subG = graph_solver.subG
        self.sources = list(graph_solver.source_dict.keys())
        self.targets = list(graph_solver.target_dict.keys())
        self.source_dict = graph_solver.source_dict
        self.target_dict = graph_solver.target_dict
        self.runinfo_df = graph_solver.runinfo_df
        self.study_id = graph_solver.study_id

    def read_runinfo_from_tsv(self, filename="./results/runinfo.tsv"):
        if os.path.exists(filename):
            self.runinfo_df = pd.read_csv(filename, sep='\t')
        else:
            print(f"Warning: {filename} does not exist!")


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



    def visualize_size_thresholds(self):
        """
        Visualize the number of nodes, edges, and % connected targets over PageRank thresholds.

        :param title: Title for the visualization.
        :param is_sign_consistent: If True, only visualize sign consistent paths.
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



    def visualize_threshold_elbowplot(self):
        """
        Visualize the threshold elbow plot of the graph. In other words, 
        the elbow indicates the best ratio size/number of connected targets.

        :param title: Title for the visualization.
        :param is_sign_consistent: If True, only visualize sign consistent paths.
        """
        thresholds_allpaths = self.runinfo_df['threshold'][self.runinfo_df['label'].str.endswith('allpaths__sc__shortest')].tolist()
        num_edges_allpaths = self.runinfo_df['num_edges'][self.runinfo_df['label'].str.endswith('allpaths__sc__shortest')].tolist()
        num_targets_allpaths = self.runinfo_df['targets_connected'][self.runinfo_df['label'].str.endswith('allpaths__sc__shortest')].tolist()
        num_edges_shortestpaths = self.runinfo_df['num_edges'][self.runinfo_df['label'].str.endswith('shortest__sc')].tolist()
        num_targets_shortestpaths = self.runinfo_df['targets_connected'][self.runinfo_df['label'].str.endswith('shortest__sc')].tolist()

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
        missing_targets_allpaths = [100-t for t in num_targets_allpaths]
        missing_targets_shortestpaths = [100-t for t in num_targets_shortestpaths]
        ax1.plot(num_edges_allpaths, missing_targets_allpaths, '-o')
        ax2.plot(num_edges_shortestpaths, [100-t for t in num_targets_shortestpaths], '-o')
        for t in thresholds_allpaths:
            ax1.text(num_edges_allpaths[thresholds_allpaths.index(t)], missing_targets_allpaths[thresholds_allpaths.index(t)], str(t))
            ax2.text(num_edges_shortestpaths[thresholds_allpaths.index(t)], missing_targets_shortestpaths[thresholds_allpaths.index(t)], str(t))
        ax1.set_ylim(0, 100)
        ax2.set_ylim(0, 100)

        ax1.set_title('Number of Edges vs missing Targets, Elbow Plot, All paths')
        ax2.set_title('Number of Edges vs missing Targets, Elbow Plot, Shortest paths')
        ax1.set_ylabel('% missing Targets')
        ax1.set_xlabel('Number of Edges')
        ax2.set_ylabel('% missing Targets')
        ax2.set_xlabel('Number of Edges')
        fig.show()

    


    def visualize_comptime(self):
        """
        Visualize the computation time of the graph.

        :param is_sign_consistent: If True, only visualize sign consistent paths.
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
        


    def visualize_degrees(self, selected_thresholds):
        """
        Visualize the degree distribution of the graph for selected thresholds.

        :param selected_thresholds: List of thresholds to visualize.
        """
        formatted_data = []
        for threshold in selected_thresholds:
            subset = self.runinfo_df[self.runinfo_df['label'] == threshold]
            if not subset.empty:
                degrees = subset['degrees'].iloc[0]
                for degree in degrees:
                    formatted_data.append([threshold, degree])

        formatted_df = pd.DataFrame(formatted_data, columns=['Threshold', 'Degree'])

        plt.figure(figsize=(10, 10))
        ax = pt.RainCloud(data = formatted_df, x = 'Threshold', y = 'Degree', palette = "Set2",
                        orient = 'v', width_box= .1, width_viol=1)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
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

        expected_files = [f'{self.study_id}__{threshold}.sif' for threshold in selected_thresholds]
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
        # intersection_series.index.names = combination

        plot(intersection_series, sort_by='cardinality', show_counts=True, show_percentages=True)
        plt.show()