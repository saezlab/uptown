
import networkx as nx
from shortest_path import SignedPathSolver, GraphSolverBase, GraphVisualizer
import time

class PageRankSolver(SignedPathSolver):
    
    def __init__(self, G):
        super().__init__(G)
        self.pagerank_values = {}
        self.is_reversed = False
        self.threshold = 0.001
    
    def reverse_graph(self):
        self.G = self.G.reverse()
        self.is_reversed = not self.is_reversed 

    def pagerank_solver(self, alpha=0.85, max_iter=100, tol=1.0e-6, nstart=None, weight='weight', dangling=None, personalize_for="source"):
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
        
        # Store the pagerank values for use in other methods
        self.pagerank_values[personalize_for] = pagerank
        
        # Set the pagerank value as an attribute for each node
        for node, pr_value in pagerank.items():
            attribute_name = 'pagerank_from_targets' if personalize_for == "target" else 'pagerank_from_sources'
            self.G.nodes[node][attribute_name] = pr_value
        
        if personalize_for == "target" and self.is_reversed:
            self.reverse_graph()

    def compute_overlap(self):
        nodes_above_threshold_from_sources = {node for node, data in self.G.nodes(data=True) if data.get('pagerank_from_sources') > self.threshold}
        nodes_above_threshold_from_targets = {node for node, data in self.G.nodes(data=True) if data.get('pagerank_from_targets') > self.threshold}
        overlap = nodes_above_threshold_from_sources.intersection(nodes_above_threshold_from_targets)
        return nodes_above_threshold_from_sources, nodes_above_threshold_from_targets, overlap
    
    def compute_filtered_shortest_paths(self, label):
        # Get nodes that exceed the threshold and overlapping nodes
        nodes_from_sources, nodes_from_targets, overlapping_nodes = self.compute_overlap()
        nodes_to_include = nodes_from_sources.union(nodes_from_targets)
        
        # Create a subgraph with only the relevant nodes
        subG = self.G.subgraph(nodes_to_include)
        
        # Store the current graph, replace with the subgraph, compute paths and then restore the original graph
        original_graph = self.G
        self.G = subG
        
        # Use the shortest_paths method from SignedPathSolver
        shortest_path, runinfo = self.shortest_paths(label = label)

        # Restore the original graph
        self.G = original_graph

        return self.shortest_paths_res, runinfo

    def visualize_pagerank(self, title="PageRankGraph", export_sif=True, is_sign_consistent=True):
        if is_sign_consistent and len(self.shortest_sc_paths_res) == 0:
            print('There were no sign consistent paths for the given perturbations and downstream effects.')
            return
        
        paths = self.shortest_sc_paths_res if is_sign_consistent else self.shortest_paths_res

        visualizer = GraphVisualizer(self)

        visualizer.visualize_graph(paths, title=title, is_sign_consistent=is_sign_consistent)

        self.to_SIFfile(paths, title=title + ".sif") if export_sif else None

    def compare_thresholds(self):
        threshold = 0.01
        targets_from_paths = []

        # get non-pagerank solvers' info
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
        targets_from_paths = {path[-1] for path in self.shortest_sc_paths_res if len(path) > 1}
        self.runinfo_dict['shortest_sc_path'].update({'targets_connected': len(targets_from_paths)})
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

            self.runinfo_dict.update(runinfo)
            self.runinfo_dict[threshold_label].update({'computation_time': computation_time})
            self.runinfo_dict[threshold_label].update({'targets_connected': len(targets_from_paths)})

            self.visualize_pagerank(threshold_label, is_sign_consistent=True, export_sif=True)
            
            print(f"Threshold computed: {threshold}")

            threshold = round(threshold - 0.001, 3)
        
    def visualize_qcplots(self, title="SizeThresholds", is_sign_consistent=True):

        visualizer = GraphVisualizer(self)

        visualizer.visualize_size_thresholds(title=title, is_sign_consistent=is_sign_consistent)
        visualizer.visualize_threshold_elbowplot(title=title, is_sign_consistent=is_sign_consistent)
        visualizer.visualize_comptime(is_sign_consistent=is_sign_consistent)

    def visualize_degrees(self, selected_thresholds, is_sign_consistent=True):
        visualizer = GraphVisualizer(self)
        visualizer.visualize_degrees(selected_thresholds, is_sign_consistent=is_sign_consistent)