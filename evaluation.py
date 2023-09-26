import pandas as pd
import os
import networkx as nx

class Eval:
    def __init__(self, dirpath):
        self.dirpath = dirpath
        self.graphs, self.metadata_df = self.get_data()

    def get_datafiles(self):
        files_list = os.listdir(self.dirpath)
        sif_files = [os.path.join(self.dirpath, f) for f in files_list if f.endswith('.sif')]
        return sif_files

    def get_data(self):
        files_list = self.get_datafiles()
        metadata_list = []
        graphs = {}

        for filename in files_list:
            G = nx.read_edgelist(filename, create_using=nx.DiGraph())
            graphs[filename] = G

            parts = filename.split('.csv')[0].split('__')
            metadata_list.append([filename] + parts)

        max_length = max(map(len, metadata_list))
        for item in metadata_list:
            while len(item) < max_length:
                item.append(None)

        metadata_df = pd.DataFrame(metadata_list, columns= ["Original_Filename"] + [f'Part_{i+1}' for i in range(max_length)])

        return graphs, metadata_df
    



