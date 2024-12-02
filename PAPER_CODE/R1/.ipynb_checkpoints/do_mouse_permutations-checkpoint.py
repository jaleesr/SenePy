import scanpy as sc
import numpy as np
import pandas as pd
import concurrent.futures
import os

import senepy as sp


# Shuffle each gene independently
def shuffle_columns(X):
    for i in range(X.shape[1]):
        np.random.shuffle(X[:, i])
    return X

def get_stats(test, adata):
    ages = ['1m', '3m', '18m', '21m', '24m', '30m']
    age_grouped = test.groupby('age', observed = True)
    
    # Create a lookup dictionary for counts and gene expression by age
    counts = age_grouped['Count'].first().to_dict()
    gene_expression = age_grouped.first().to_dict()
    
    
    out = []
    for gene in adata.var_names:
        if gene == 'age' or gene == 'Count':
            continue
        
        age_data = {}
        for age in ages:
            age_data[f'c{age}'] = counts.get(age, 0)
            age_data[f'f{age}'] = gene_expression.get(gene, {}).get(age, 0)
        
        # Compute young_start and old_end
        if age_data['c3m'] > 200:
            young_start = age_data['f3m']
        else:
            young_start = (age_data['c3m'] * age_data['f3m'] + age_data['c1m'] * age_data['f1m']) / (age_data['c3m'] + age_data['c1m'])


        #need fix for mammary because mammary is the only one without 24m and 30m, using 21m instead
        if age_data['c30m'] > 200:
            old_end = age_data['f30m']
        elif test.iloc[0].tissue2 != 'Mammary_Gland':
            old_end = (age_data['c30m'] * age_data['f30m'] + age_data['c24m'] * age_data['f24m']) / (age_data['c30m'] + age_data['c24m'])
        else:
            old_end = (age_data['c21m'] * age_data['f21m']) / (age_data['c21m'])
    
        # Calculate old_gain and old_ratio
        old_gain = old_end - young_start
        old_ratio = old_end / young_start if young_start != 0 else 0
        
        # Append results to output
        out.append([gene, young_start, old_end, old_gain, old_ratio])
    
    # Convert the output to DataFrame
    out_df = pd.DataFrame(out, columns=['gene', 'start', 'end', 'gain', 'ratio'])

    return out_df


def do_permutation(tissue, cell, perms):

    f = f'temp_adata_subs/{tissue}_{cell}.h5ad'
    f = f.replace(' ', '_')
    
    print(f"Starting permutations for tissue: {tissue}, cell: {cell}")

    for ii in range(perms):
        if os.path.exists(f'permutation_output/{tissue}.{cell}.{ii}.csv'):
            continue

        cdata = sc.read_h5ad(f)
        
        cdata.X = cdata.X.toarray()
        
        sc.pp.filter_genes(cdata, min_cells=1)
        
        cdata.X = shuffle_columns(cdata.X)
        
        cell_counts = cdata.obs.groupby(['age', 'tissue2', 'cell_type_2'], observed=True)\
        .size().reset_index().rename(columns={0:'Count'})
        
        cell_counts = cell_counts[cell_counts.Count > 0]
        
        out = []
        for age in cell_counts.age.values:
            age_sub = cdata[cdata.obs.age == age].X.toarray()
            
            col_out = []
            for gene in cdata.var_names:
                i = np.where(cdata.var_names == gene)[0][0]
                gene_vals = age_sub[:, i]
                gene_frac = len(gene_vals[gene_vals > 0])/len(gene_vals)
                if gene_frac == 0 and age in ['1m', '3m']: #only impute young
                    gene_frac = 1/len(gene_vals) #impute 1/num cells if value is 0
                
                col_out.append(gene_frac)
            out.append(col_out)
            
        out = np.array(out)
    
        cell_counts = pd.concat([cell_counts, pd.DataFrame(out, columns = list(cdata.var_names))], axis = 1)
    
        
        the_stats = get_stats(cell_counts, cdata)
    
        the_stats.to_csv(f'permutation_output/{tissue}.{cell}.{ii}.csv')

    print(f"Done permutations for tissue: {tissue}, cell: {cell}")




def main():
    N_perms = 1000

    # ####FOR SENEPY SIGNATURES, the main reason for the script
    # hubs = sp.load_hubs(species='Mouse', sig_type='cell_type')
    # input_list = [list(x) + [N_perms] for x in list(hubs.hubs)]

    ####this is run and above two lines commented out only when running for additional perms for fig 1
    amp = pd.read_csv('additional_mouse_perms.csv', index_col=0)
    input_list = [list(x) + [N_perms] for x in zip(amp.tissue, amp.cell_type_2)]


    
    with concurrent.futures.ProcessPoolExecutor(max_workers=40) as executor:
        futures = [executor.submit(do_permutation, *x) for x in input_list]
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"Error occurred: {e}")

if __name__ == "__main__":
    main()
















