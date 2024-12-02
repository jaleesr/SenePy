import scanpy as sc
import numpy as np
import pandas as pd
import concurrent.futures
import os
from scipy import stats

import senepy as sp


# Shuffle each gene independently
def shuffle_columns(X):
    for i in range(X.shape[1]):
        np.random.shuffle(X[:, i])
    return X

def get_stats(x):
    ages = np.array([8,18,28,38,48,58,68,78,88])
    mask = x > -1
    
    ages = ages[mask]
    x = x[mask]
    
    lin = stats.linregress(ages, x)
    
    if 18 in ages or 28 in ages or 8 in ages:
        starting = x.iloc[0]
    else:
        starting = lin.slope*18 + lin.intercept
        
    if starting < 0:
        starting = 0
        
    ending = x.iloc[-1]

    return lin.slope, x.iloc[0], x.max(), lin.pvalue, len(x), starting, ending


def do_permutation(tissue, cell, perms):

    f = f'human_adatas/{tissue}.{cell}.h5ad'
    
    print(f"Starting permutations for tissue: {tissue}, cell: {cell}")

    for ii in range(perms):
        if os.path.exists(f'permutation_output_human/{tissue}.{cell}.{ii}.csv'):
            continue

        cdata = sc.read_h5ad(f)
        
        cdata.X = cdata.X.toarray()
        
        sc.pp.filter_genes(cdata, min_cells=1)
        svars = cdata.var_names
        
        cdata.X = shuffle_columns(cdata.X)


        temp = []
        for age in [8,18,28,38,48,58,68,78,88]:
            sdata = cdata[cdata.obs.bin_age == age].X
            for gene in svars:
                temp2 = [age, gene]
                if len(sdata) > 100:
                    i = np.where(svars == gene)[0][0]
                    a = sdata[:,i]
                    v = a[a > 0].shape[0] / a.shape[0]
                else:
                    v = -1
                temp2 += [v]
                temp.append(temp2)
        
        tc_df = pd.DataFrame(temp, columns = ['age_bin', 'gene', 'value'])
        tc_df = tc_df.pivot(values = 'value', index = ['gene'], columns = 'age_bin').reset_index()

        
        tc_df['slope'], tc_df['starting'], tc_df['Max'], tc_df['lin_p'], tc_df['N'], tc_df['starting'], tc_df['ending'] =\
                zip(*tc_df[[8,18,28,38,48,58,68,78,88]].apply(get_stats, axis = 1))

    
        tc_df.to_csv(f'permutation_output_human/{tissue}.{cell}.{ii}.csv')

    print(f"Done permutations for tissue: {tissue}, cell: {cell}")




def main():
    N_perms = 1000
    # adata_path = '../../data/tms/tms-scVI-raw-data_BDATA.h5ad'
    # #adata = sc.read_h5ad(adata_path)
    hubs = sp.load_hubs(species='Human', sig_type='cell_type')

    input_list = [list(x) + [N_perms] for x in list(hubs.hubs)]

    with concurrent.futures.ProcessPoolExecutor(max_workers=40) as executor:
        futures = [executor.submit(do_permutation, *x) for x in input_list]
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"Error occurred: {e}")

if __name__ == "__main__":
    main()














