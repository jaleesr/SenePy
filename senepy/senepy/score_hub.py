import numpy as np
import pandas as pd
from tqdm import tqdm

def score_hub(adata,
              hub,
              n_bins = 25,
              ctrl_size = 50,
              binarize = True,
              importance = True,
              translator = None,
              verbose = True):
    
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
    var_names = adata.var_names
    
    
    
    genes_not_in_adata = [x[0] for x in hub if x[0] not in var_names]
    
    if len(genes_not_in_adata) > 0 and verbose:
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
        
        
        if verbose:
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
    
    
    
    
    
def score_all_cells(adata,
              hub,
              identifiers = None,
              n_bins = 25, 
              ctrl_size = 50,
              binarize = True,
              importance = True,
              translator = None):
    
    
    
    '''
    Scores all cell types. Takes anndata, a hub list, and .obs identifier columns. Breaks
    anndata into chunks based on identifiers then scores them individiaully. Returns a
    score list for all cells.
    
    adata: anndata object
    
    hub: list
        list of tuples: [(gene, importance), (gene, importance), ... ]
        
    identifiers: list
        list of columns in .obs to split the data on. Works with multiple:
        e.g., ['tissue', 'cell_type'] would seperate lung fibroblasts from heart
        fibroblasts
        
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
        
    
    '''
    
    if identifiers is None and isinstance(identifiers, list):
        raise ValueError("Provide cell_type identifiers (column names) from the .obs dataframe as a list")
        
    #cdata = adata.copy()
    
    var_names = adata.var_names
    
    genes_not_in_adata = [x[0] for x in hub if x[0] not in var_names]
   
            
    unique_combos = adata.obs[identifiers].drop_duplicates()
            
    score_d = {}
    to_be_verbose = True
    for idents in tqdm(unique_combos.values):
        
        #subset adata based on identifiers. For each id, keep only those
        c_sub = adata
        for x, ident in enumerate(identifiers):
            c_sub = c_sub[c_sub.obs[ident] == idents[x]]
        
        
        
        res = score_hub(c_sub, hub, n_bins = n_bins, ctrl_size = ctrl_size,
                  binarize = binarize, importance = importance,
              translator = translator, verbose = to_be_verbose)
            
        score_d = score_d | dict(zip(c_sub.obs.index, res))
        
        to_be_verbose = False
        
        
    
    return adata.obs.index.map(score_d).tolist()
    
    


def score_hub_scdrs(adata,
                    hub,
                    stats_key = 'scdrs_stats',
                    translator = None,
                    verbose = True,
                    citation = False,
                    **kwargs):   

    '''
    Wrapper for scdrs that scores all cells based on an input senePy hub. Returns a list of
    scores and updates the input anndata with a statistics dataframe in .uns.
    
    Retains most functionality from scdrs
    
    adata: anndata object
        .X should be normalized and transformed
    
    hub: list
        list of tuples: [(gene, importance), (gene, importance), ... ]
        
    translator: translator object
        optional translator object which will map the hub genes to have better overlap
        with your study genes. See senepy.translator
        
    stats_key: str
        Name of the key to add in anndata.uns for the resulting stats dataframe that will
        be added to your anndata object.
    
     Refer to `scdrs.pp.preprocess` and `scdrs.score_cell` for more details and additonal
     arguments that can be used.
    
    
    Reference:
        This is a modified approach from Zhang et al. (2022),
        Polygenic enrichment distinguishes disease associations of individual cells
        in single-cell RNA-seq data, Nature Genetics. 
        
    
    '''

    
    try:
        import scdrs
    except ImportError:
        raise ImportError("scDRS is required to use the score_hub_scdrs. "
                              "Please install it by running 'pip install scdrs==1.0.2'.")

    if not citation:
        print('Please cite scDRS if you use this function in your work: doi.org/10.1038/s41588-022-01167-z')
        print('To quiet this message pass citation = True')



    #translation of var names, should maybe make this into its own function since I have now done it three times...
    var_names = adata.var_names
    genes_not_in_adata = [x[0] for x in hub if x[0] not in var_names]
    
    if len(genes_not_in_adata) > 0 and verbose:
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
        
        
        if verbose:
            print(f'{len(hub)-len(genes_not_in_adata)}/{len(hub)}({round((len(hub)-len(genes_not_in_adata))/len(hub)*100,2)}%) genes present in data after translation')
            print('Still not present:', genes_not_in_adata)


    
    #scDRS scoring...
    defaults = {'cov': None, 'adj_prop': None, 'n_mean_bin': 20, 'n_var_bin': 20, 'n_chunk': None, 
               'ctrl_match_key': "mean_var", 'n_ctrl': 1000, 'n_genebin': 200, 'weight_opt': 'vs'}
    
    options = {**defaults, **kwargs}

    if adata.X.min() < 0 or adata.X.max() > 30:
        raise ValueError("The AnnData object is not properly normalized and trasnformed. Please provied normalized-transfored data")


    pdata = scdrs.pp.preprocess(adata, cov = options['cov'],
                        adj_prop = options['adj_prop'],
                        n_mean_bin = options['n_mean_bin'],
                        n_var_bin = options['n_var_bin'],
                        n_chunk = options['n_chunk'],
                        copy = True)


    genes = [x[0] for x in hub]
    weights = [x[1] for x in hub]


    res = scdrs.score_cell(
                pdata,
                genes,
                gene_weight = weights,
                ctrl_match_key = options['ctrl_match_key'],
                n_ctrl = options['n_ctrl'],
                n_genebin = options['n_genebin'],
                weight_opt = options['weight_opt'],
                copy=False,
                return_ctrl_raw_score=False,
                return_ctrl_norm_score=False,
                random_seed=0,

            )


    adata.uns[stats_key] = res

    return res.norm_score.tolist()
    











    
    
    
    
    
    
    