import numpy as np
import pandas as pd


def score_hub(adata,
              hub,
              n_bins = 25,
              ctrl_size = 50,
              binarize = True,
              importance = True,
              translator = None):
    
    '''
    Takes anndata and a hub list and returns a score for each cell
    
    adata: anndata object
    
    hub: list
        list of tuples: [(gene, importance), (gene, importance), ... ]
        
    n_bins: int
        number of expression bins to split the genes by. The expression for each gene is 
        averaged then they are ordered and split into n bins. The background genes are
        selected from the same bins the input genes fall in.
        
    ctrl_size: int
        number of background genes for each input gene
        
    binarize: bool
        whether to binarize the data after backround gene selection. 1 if the gene is
        expressed and 0 if it is not.
        
    importance: bool
        whether to use the importance associated with each gene. Importance may be less
        relevant if the cells you are testing are much different from the hub cells.
        
    translator: translator object
        optional translator object which will map the hub genes to have better overlap
        with your study genes. See senepy.translator
        
        
        
    Reference:
        This is a modified approach from Satija et al. (2015),
        Spatial reconstruction of single-cell gene expression data, Nature Biotechnology. 
    
    
    '''
    
    
    cdata = adata.copy()
    
    np.random.seed(0)
    var_names = cdata.var_names
    
    
    
    genes_not_in_adata = [x[0] for x in hub if x[0] not in var_names]
    
    if len(genes_not_in_adata) > 0:
        print(f'{len(hub)-len(genes_not_in_adata)}/{len(hub)}({round((len(hub)-len(genes_not_in_adata))/len(hub)*100,2)}%) genes present in data')
        if translator is None:
            print('###################')
            print('Not present:', genes_not_in_adata)
            print('###################')
            print('passing a translator may improve overlap')
            
        
    #if a translator is passed, convert the hub genes based on the translator.mapper
    if translator is not None:
        hub_cpy = []
        for gene, y in hub:
            try:
                hub_cpy.append((translator.mapper[gene], y))
            except:
                hub_cpy.append((gene, y))
    
        hub = hub_cpy
    
        genes_not_in_adata = [x[0] for x in hub if x[0] not in var_names]
        print(f'{len(hub)-len(genes_not_in_adata)}/{len(hub)}({round((len(hub)-len(genes_not_in_adata))/len(hub)*100,2)}%) genes present in data after translation')
        print('Still not present:', genes_not_in_adata)
    
    
    
    
    hub = [x for x in hub if x[0] in var_names]
    
    present_genes = [x[0] for x in hub]
    
    if len(present_genes) == 0:
        raise ValueError("No genes matched your dataset. Try using a translator")

    
    #densify to make amplification faster
    try:
        cdata.X = cdata.X.toarray() #faster but more memory
    except:
        pass
    
    
    
    
    ######## time to select background genes ########
    
    #get mean expression for each gene
    gene_exp_avg = pd.Series(np.nanmean(cdata.X, axis=0), index = cdata.var_names)
    gene_exp_avg = gene_exp_avg[np.isfinite(gene_exp_avg)] #sometimes data missing?
    
    
#     #if i wanted to rank by dropout instead
#     nash.X.getnnz(axis = 0) #get num of nonzero per gene, FOR NONSPARSE
#     np.count_nonzero(nash.X.toarray(), axis=0) #for dense
    
    
    # num of genes -1 / n_bins
    n_items = int(np.round(len(gene_exp_avg) / (n_bins - 1)))
    
    
    
    #first gives a numeric rank from min to the gene averages, then floor divide by previous n_items
    gene_ranks = gene_exp_avg.rank(method='min') // n_items
    
    
    
    control_genes = set()
    # now pick `ctrl_size` genes from every gene rank
    for rank in np.unique(gene_ranks.loc[present_genes]): 
        r_genes = np.array(gene_ranks[gene_ranks == rank].index) #genes with equal rank values to input gene
        #len r_genes should equal n_tems for every rank
        np.random.shuffle(r_genes)
        control_genes.update(set(r_genes[:ctrl_size])) #takes the first ctrl_size num of shuffled genes
        
    ######## done background selection ########
    
    
    
    ### binarize
    if binarize:
        cdata.X[cdata.X > 0] = 1
        
    
    
    #### amplify expression matrix by gene importance
    if importance:
        for gene, importance in hub:
            i = np.where(cdata.var_names == gene)[0][0]
            cdata.X[:,i] = cdata.X[:,i] * importance
        
        
        
        
     #### time to get scores ####

    #mean for input genes for each cell
    X_list = cdata[:, present_genes].X #cells X input genes matrix
    X_list = np.nanmean(X_list, axis=1, dtype='float64')

    #X_control is mean of control genes for each cell
    control_genes = list(control_genes - set(present_genes))
    X_control = cdata[:, control_genes].X #cells X control genes
    X_control = np.nanmean(X_control, axis=1, dtype='float64')
        
    

    score = X_list - X_control

    
    return score.tolist()
    
    
    