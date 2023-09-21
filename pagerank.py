
import networkx as nx
from shortest_path import BaselineSolver, GraphSolverCore, GraphVisualizer
import time


class PageRankSolver(BaselineSolver):
    """
    Solver for computing PageRank values and related operations on a graph.
    A child class of the parent class BaselineSolver.
    """
    
    def __init__(self, G, study_id = None):
        """
        Initialize the PageRank solver with a given graph.

        :param G: A networkx graph object.
        """
        super().__init__(G, study_id)
        self.pagerank_values = {}
        self.is_reversed = False
        self.threshold = 0.001
    


    def reverse_graph(self):
        """
        Reverse the direction of all edges in the graph.
        """
        self.G = self.G.reverse()
        self.is_reversed = not self.is_reversed 



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
        return nodes_above_threshold_from_sources, nodes_above_threshold_from_targets, overlap
    


    def compute_filtered_shortest_paths(self, label):
        """
        Compute the shortest paths in the graph after filtering nodes based on PageRank threshold.

        :param label: A label for the current run.
        :return: A tuple containing the shortest paths and run information.
        """
        # Get nodes that exceed the threshold and overlapping nodes
        nodes_from_sources, nodes_from_targets, overlapping_nodes = self.compute_overlap()
        nodes_to_include = nodes_from_sources.union(nodes_from_targets)
        
        # Create a subgraph with only the relevant nodes
        subG = self.G.subgraph(nodes_to_include)
        
        # Store the current graph, replace with the subgraph, compute paths and then restore the original graph
        original_graph = self.G
        self.G = subG

        shortest_path, runinfo = self.shortest_paths(label = label)

        # Restore the original graph
        self.G = original_graph

        return self.all_paths_res, runinfo
    
    def compute_filtered_all_paths(self, label):
        """
        Compute the shortest paths in the graph after filtering nodes based on PageRank threshold.

        :param label: A label for the current run.
        :return: A tuple containing the shortest paths and run information.
        """
        # Get nodes that exceed the threshold and overlapping nodes
        nodes_from_sources, nodes_from_targets, overlapping_nodes = self.compute_overlap()
        nodes_to_include = nodes_from_sources.union(nodes_from_targets)
        
        # Create a subgraph with only the relevant nodes
        subG = self.G.subgraph(nodes_to_include)
        
        # Store the current graph, replace with the subgraph, compute paths and then restore the original graph
        original_graph = self.G
        self.G = subG

        shortest_path, runinfo = self.all_paths(label = label)

        # Restore the original graph
        self.G = original_graph

        return self.shortest_paths_res, runinfo



    def visualize_pagerank(self, title="PageRankGraph", export_sif=True, is_sign_consistent=True):
        """
        Wrapper around the visualize_graph method from GraphVisualizer class. Provides a graph view of the network.

        :param title: Title for the visualization.
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



    def compare_thresholds(self, threshold = 0.01):
        """
        Computes shortest paths and sign consistent paths for different PageRank thresholds.
        """
        targets_from_paths = []

        start_time = time.time()
        shortest_paths, runinfo = self.shortest_paths(label = 'shortest_path')
        end_time = time.time()
        computation_time = end_time - start_time
        self.runinfo_dict.update(runinfo)
        self.runinfo_dict['shortest_path'].update({'computation_time': computation_time})
        targets_from_paths = {path[-1] for path in self.shortest_paths_res if len(path) > 1}
        self.runinfo_dict['shortest_path'].update({'targets_connected': len(targets_from_paths)})
        self.visualize_graph("shortest_path", is_sign_consistent=False, export_sif=True)

        start_time = time.time()
        shortest_sc_paths, runinfo = self.sign_consistency_check('shortest_sc_path')
        end_time = time.time()
        computation_time = end_time - start_time
        self.runinfo_dict.update(runinfo)
        self.runinfo_dict['shortest_sc_path'].update({'computation_time': computation_time})
        sc_targets_from_paths = {path[-1] for path in self.shortest_sc_paths_res if len(path) > 1}
        self.runinfo_dict['shortest_sc_path'].update({'targets_connected': len(sc_targets_from_paths)})
        self.visualize_graph("shortest_sc_path", is_sign_consistent=True, export_sif=True)

        while len(targets_from_paths) < len(self.target_dict) and threshold >= 0:
            
            self.threshold = threshold
            threshold_label = f"pagerank_{threshold}"

            start_time = time.time()
            paths, runinfo = self.compute_filtered_shortest_paths(label = threshold_label)
            end_time = time.time()

            targets_from_paths = {path[-1] for path in self.shortest_paths_res if len(path) > 1}

            self.runinfo_dict.update(runinfo)

            self.runinfo_dict[threshold_label].update({'computation_time': computation_time})
            self.runinfo_dict[threshold_label].update({'targets_connected': len(targets_from_paths)})

            self.visualize_pagerank(threshold_label, is_sign_consistent=False, export_sif=True)


            threshold_label = f"pagerank_sc_{threshold}"

            start_time = time.time()
            paths, runinfo = self.sign_consistency_check(label = threshold_label)
            end_time = time.time()

            computation_time = end_time - start_time

            sc_targets_from_paths = {path[-1] for path in self.shortest_sc_paths_res if len(path) > 1}

            self.runinfo_dict.update(runinfo)
            self.runinfo_dict[threshold_label].update({'computation_time': computation_time})
            self.runinfo_dict[threshold_label].update({'targets_connected': len(sc_targets_from_paths)})

            self.visualize_pagerank(threshold_label, is_sign_consistent=True, export_sif=True)
            
            print(f"Threshold computed: {threshold}")

            threshold = round(threshold - 0.001, 3)
        


    def visualize_qcplots(self, title="SizeThresholds", is_sign_consistent=True):
        """
        Visualize quality control plots for the graph. Wrapper around the 
        visualize_size_thresholds, visualize_threshold_elbowplot and visualize_comptime 
        methods from GraphVisualizer class.

        :param title: Title for the visualization.
        :param is_sign_consistent: If True, only visualize sign consistent paths.
        """
        visualizer = GraphVisualizer(self)

        visualizer.visualize_size_thresholds(title=title, is_sign_consistent=is_sign_consistent)
        visualizer.visualize_threshold_elbowplot(title=title, is_sign_consistent=is_sign_consistent)
        visualizer.visualize_comptime(is_sign_consistent=is_sign_consistent)



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