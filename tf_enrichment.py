import decoupler as dc
import pandas as pd

panacea_de = pd.read_csv('panacea_de.tsv', sep='\t')
panacea_de.drop(columns=['logFC', 'padj'], inplace=True)
# join cell line and treatment columns
panacea_de['bio_context'] = panacea_de['cell_line'].str.cat(panacea_de['treatment'], sep='_')
panacea_de.drop(columns=['cell_line', 'treatment'], inplace=True)
panacea_de = panacea_de.pivot_table(index='gene_symbol', columns='bio_context', values='stat').transpose()
panacea_de.fillna(0, inplace=True)
mat = panacea_de.astype(float)
#convert to matrix


# decoupler process
net = pd.read_csv('collectri.tsv', sep='\t')
net
net.drop(columns = 'PMID', inplace=True)


results = dc.dense_run(dc.run_ulm, mat, net)
t_values, p_values = results
t_values.to_csv('tf_activity_results.tsv', sep='\t')