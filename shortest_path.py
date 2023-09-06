# shortest_path.py

import pandas as pd
import networkx as nx
import pygraphviz as pgv
import numpy as np

class GraphSolverBase:
    def __init__(self, G):
        self.G = G

class SignedPathSolver(GraphSolverBase):

    def __init__(self, G):
        super().__init__(G)
        self.source_dict = {}
        self.target_dict = {}
        self.shortest_paths_res = []
        self.shortest_sc_paths_res = []
        self.connected_targets = {}
        self.connected_sc_targets = {}

    def shortest_paths(self, verbose = False):
        self.shortest_paths_res = []
        self.connected_targets = {}
        sources = list(self.source_dict.keys())
        targets = list(self.target_dict.keys())
        for source_node in sources:
            # Initialize an empty list for the source_node if it doesn't exist
            if source_node not in self.connected_targets:
                self.connected_targets[source_node] = []

            for target_node in targets:
                try:
                    self.shortest_paths_res.extend([p for p in nx.all_shortest_paths(self.G, source=source_node, target=target_node, weight='weight')])
                    self.connected_targets[source_node].append(target_node)
                except nx.NetworkXNoPath as e:
                    if verbose:
                        print(f"Warning: {e}")

    def sign_consistency_check(self):
        self.shortest_sc_paths_res = []
        self.connected_sc_targets = {}
        for path in self.shortest_paths_res:
            product_sign = 1  # initialize the product of signs to 1
            source = path[0]
            target = path[-1]
            source_sign = self.source_dict[source]
            target_sign = self.target_dict[target]
        
            if source not in self.connected_sc_targets:
                self.connected_sc_targets[source] = []

            # calculate the product of the signs of the edges
            for i in range(len(path) - 1):
                edge_sign = self.G.get_edge_data(path[i], path[i + 1])['sign']
                product_sign *= edge_sign
            
            # check if the product of the signs matches the expectation
            if product_sign == source_sign * target_sign:
                self.shortest_sc_paths_res.append(path)  # append the consistent path to sscps list
                self.connected_sc_targets[source].append(target) # append the target to the list of connected targets
            self.connected_sc_targets[source] = list(set(self.connected_sc_targets[source]))

    
    def to_SIFfile(self, paths, title):
        sif_tuples = []
        for path in paths:
            for i in range(len(path) - 1):
                interaction_type = 'P' if self.G[path[i]][path[i+1]]['sign'] > 0 else 'N'
                sif_tuples.append((path[i], interaction_type, path[i+1]))
        
        sif_df = pd.DataFrame(sif_tuples, columns=['source', 'interaction', 'target']).drop_duplicates()
        sif_df.to_csv(title, sep='\t', index=None)
        return(sif_df)

    def visualize_graph(self, paths, title="Graph", export=False, is_sign_consistent=True):
        if is_sign_consistent and len(self.shortest_sc_paths_res) == 0:
            print('There were no sign consistent paths for the given perturbations and downstream effects.')
            return
        
        paths = self.shortest_sc_paths_res if is_sign_consistent else self.shortest_paths_res

        visualizer = GraphVisualizer(self)

        title = title + "_sign_consistent" if is_sign_consistent else title

        visualizer.visualize_graph(paths, title=title, is_sign_consistent=is_sign_consistent)  



class GraphVisualizer:
    def __init__(self, graph_solver):
        self.G = graph_solver.G
        self.sources = list(graph_solver.source_dict.keys())
        self.targets = list(graph_solver.target_dict.keys())
        self.source_dict = graph_solver.source_dict
        self.target_dict = graph_solver.target_dict

    def visualize_graph(self, paths, title="Graph", export=False, is_sign_consistent=True):
        

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

        # Convert networkx graph to pygraphviz graph
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

        # Add styles to nodes
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
        
        # Add styles to edges
        for edge in A.edges():
            u, v = edge
            edge_data = V.get_edge_data(u, v)
            edge_color = 'green' if edge_data['sign'] == 1 else 'red'
            edge.attr['color'] = edge_color
        
        # Display the graph
        A.layout(prog='dot')
        file_path = f'{title}.pdf'
        A.draw(file_path, prog='dot')

        # Export file
        self.to_SIFfile(paths, title + ".sif") if export else None