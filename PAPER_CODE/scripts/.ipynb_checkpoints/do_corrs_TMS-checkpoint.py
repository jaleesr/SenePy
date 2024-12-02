import scanpy as sc
import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm
import itertools
import concurrent.futures
import sys
import time
import os
from numba import jit

#note that Pearsons'R simplifies and aligns to Phi

@jit(nopython=True)
def fast_perm(g1, g2):
    the_corr = np.corrcoef(g1, g2)[:,0][1]
    out = [1]
    for x in range(500):
        np.random.shuffle(g2)
        rnd_corr = np.corrcoef(g1, g2)[:,0][1]
        out.append(rnd_corr)   
        
    out = np.array(out[1:])
    p = len(out[out >= the_corr] )/500
    
    return [the_corr, p, np.quantile(out, .99)]


fast_perm(np.array([0,1,1,1,1,0]), np.array([0,1,1,1,1,0])) #initialize jit



def a_function(chunk):
    out = []
    for combo in chunk:
        i1 = np.where(sub.var_names == combo[0])[0][0]
        i2 = np.where(sub.var_names == combo[1])[0][0]
        a1 = sub_X[:, i1].copy()
        a2 = sub_X[:, i2].copy()
        a1[a1 > 0] = 1 #binarize
        a2[a2 > 0] = 1 #binarize
        
        perms = fast_perm(a1.copy(), a2.copy())
        

        out.append([tissue,cell,combo[0], combo[1]] + perms)
    return out
    
    
    
    
threads = int(sys.argv[2])
sub = sc.read_h5ad(sys.argv[1])
sub_X = sub.X.toarray()

tissue = sub.obs.tissue2.iloc[1]
cell = sub.obs.cell_type_2.iloc[1]


out_str = 'out_dir/' + tissue + '__' + cell + '.pickle'
out_str = out_str.replace(' ', '_')


doners = os.listdir('out_dir')
if out_str.split('/')[1] in doners:
    sys.exit(f"{out_str} already exists")



the_combos = list(itertools.combinations(sub.var_names, 2))



chunks = []
for x in range(0, len(the_combos), 200):
    chunks.append(the_combos[x:x+200])

print(f'starting {tissue} : {cell}')    
print(f'# {len(the_combos)} total combos and {len(chunks)} chunks')

start = time.time()

with concurrent.futures.ProcessPoolExecutor(max_workers = threads) as executor:
    iterator = list(tqdm(executor.map(a_function, chunks)))
    

end = time.time()
print(end-start)
    

broken = []
for item in list(iterator):
    broken += item
    
corrs = pd.DataFrame(broken, columns = ['tissue', 'cell_type_2', 'gene1', 'gene2', 'r', 'p', 'rnd_r_q99'])

    


corrs.to_pickle(out_str)
































