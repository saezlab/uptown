from matplotlib import cm
import networkx as nx
from shortest_path import SignedPathSolver, GraphSolverBase, GraphVisualizer

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
    
    def compute_filtered_shortest_paths(self, sign_consistent=True):
        # Get nodes that exceed the threshold and overlapping nodes
        nodes_from_sources, nodes_from_targets, overlapping_nodes = self.compute_overlap()
        nodes_to_include = nodes_from_sources.union(nodes_from_targets)
        
        # Create a subgraph with only the relevant nodes
        subG = self.G.subgraph(nodes_to_include)
        
        # Store the current graph, replace with the subgraph, compute paths and then restore the original graph
        original_graph = self.G
        self.G = subG
        
        # Use the shortest_paths method from SignedPathSolver
        self.shortest_paths()
        
        # Optionally filter paths based on sign consistency
        if sign_consistent:
            self.sign_consistency_check()
        
        # Restore the original graph
        self.G = original_graph
        print(len(overlapping_nodes))

        return self.shortest_sc_paths_res if sign_consistent else self.shortest_paths_res

    def visualize_pagerank(self, title="PageRankGraph", is_sign_consistent=True):
        if is_sign_consistent and len(self.shortest_sc_paths_res) == 0:
            print('There were no sign consistent paths for the given perturbations and downstream effects.')
            return
        
        title = title + "_sign_consistent" if is_sign_consistent else title
        
        paths = self.shortest_sc_paths_res if is_sign_consistent else self.shortest_paths_res

        visualizer = GraphVisualizer(self)

        visualizer.visualize_graph(paths, title=title, is_sign_consistent=is_sign_consistent)  