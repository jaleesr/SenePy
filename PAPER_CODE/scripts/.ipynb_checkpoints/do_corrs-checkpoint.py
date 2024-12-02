import scanpy as sc
import numpy as np
import pandas as pd

from scipy import stats

from tqdm import tqdm

import itertools
import concurrent.futures

import sys
import os

#note that Pearsons'R simplifies and aligns to Phi

doners = os.listdir('out_dir/')



sub = sc.read_h5ad(sys.argv[1])
tissue = sub.obs.tissue2.iloc[1]
cell = sub.obs.cell_type_2.iloc[1]


out_str = 'out_dir/' + tissue + '__' + cell + '.pickle'
out_str = out_str.replace(' ', '_')


if out_str.split('/')[1] in doners:
    sys.exit(f"{out_str} already exists")
    
    



print(f'tissue: {tissue}\t cell: {cell}')


cell_stats = pd.read_pickle('/home/ubuntu/s3_mount/senescence/paper_1/data/files/tms_corr_subs/the_stats.pickle')


the_genes = cell_stats[(cell_stats.tissue2 == tissue) &\
                    (cell_stats.cell_type_2 == cell)].gene.tolist()




# tissue, cell = aa

sub_X = sub.X.toarray()

out = []
for combo in tqdm(list(itertools.combinations(the_genes, 2))):
    i1 = np.where(sub.var_names == combo[0])[0][0]
    i2 = np.where(sub.var_names == combo[1])[0][0]
    a1 = sub_X[:, i1].copy()
    a2 = sub_X[:, i2].copy()
    a1[a1 > 0] = 1 #binarize
    a2[a2 > 0] = 1 #binarize
    v = stats.pearsonr(a1, a2)


    ps = []
    rs = []

#     if len(the_genes) < 100:
#         perms = 1000
#     elif len(the_genes) < 750:
#         perms = 500
#     else:
#         perms = 200

    for x in range(0,500): #perms
        np.random.shuffle(a1)
        np.random.shuffle(a2)
        vc = stats.pearsonr(a1, a2) #r, p
        ps.append(vc[1])
        rs.append(vc[0])

    ps = np.array(ps)
    rs = np.array(rs)

    out.append([tissue,cell,
                combo[0], combo[1], v[0], v[1],
                rs.mean(), ps.mean(),
                np.quantile(rs, .99), np.quantile(ps, .01), rs.std(), ps.std()])



test = pd.DataFrame(out, columns = ['tissue', 'cell_type_2', 'gene1', 'gene2', 'r', 'p',
                             'rnd_r_mean', 'rnd_p_mean', 'rnd_r_q99', 'rnd_p_q01', 'rnd_r_sd', 'rnd_p_sd'])



test.to_pickle(out_str)
































