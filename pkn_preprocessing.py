import omnipath as op
import numpy as np
import pandas as pd
import decoupler as dc

collectri = dc.get_collectri()
omnipath = op.interactions.AllInteractions.get('omnipath',genesymbols=True)

omnipath.rename(columns={'source':'source_uniprot', 'target':'target_uniprot'}, inplace=True)
omnipath.rename(columns={'source_genesymbol':'source', 'target_genesymbol':'target'}, inplace=True)

# get only directed interactions
omnipath = omnipath[(omnipath['consensus_direction'] == True) & (((omnipath['consensus_stimulation'] == True) & (omnipath['consensus_inhibition'] == False)) | ((omnipath['consensus_stimulation'] == False) & (omnipath['consensus_inhibition'] == True))) & (omnipath['curation_effort'] >= 2)]

# remove interactions downstream of TFs
collectri_tfs = collectri.source.unique().tolist()
omnipath = omnipath[~omnipath['source'].isin(collectri_tfs)]


# subset the collectri GRN to enrich only TFs present in the network
collectri_tfs_in_omnipath = set(omnipath['target'][omnipath['target'].isin(collectri_tfs)])
collectri_tfs_in_omnipath_df = collectri[collectri['source'].isin(collectri_tfs_in_omnipath)]
collectri_tfs_in_omnipath_df.to_csv('collectri.tsv', sep='\t', index=False)
collectri_tfs_in_omnipath_df

omnipath['data'] = np.where(omnipath['consensus_stimulation'] == True, 1, -1)

# write the resulting omnipath network in networkx format
net = omnipath[['source', 'target', 'data']]
net.rename(columns={'source':'# source'}, inplace=True)
net.to_csv('network_collectri.sif', sep = '\t', index = False)