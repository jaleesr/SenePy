import pickle
from pandas import read_pickle, DataFrame
from pkg_resources import resource_filename
from math import comb
from collections import Counter
from tqdm import tqdm
from random import choices
from numpy import array


def hypergeom_sf(x, M, n, N):
    '''
    Calculate the survival function of the hypergeometric distribution

    Used below in load_hubs.search_hubs_by_genes
    '''
    cdf = 0.0
    for k in range(0, x + 1):
        cdf += (comb(n, k) * comb(M - n, N - k)) / comb(M, N)
    sf = 1 - cdf
    return sf


class load_hubs(object):
    '''
    Load the senescence gene hubs generated in the study
    
        species: string
            'Human' or 'Mouse'
            
        sig_type: string
            'hubs' or 'cell_type': If 'hubs', signatures are split into individual hubs
                with mulitple hubs for some cell type. If 'cell_type', hubs are
                combined into one signature for individual cell types.
    
    
        self.hubs: dictionary of the hubs with genes, importances
        self.metadata: brief metadata of the hubs
            hub_num: there are multiple hubs in some individiual cell types
            size: number of genes
            n_sen: number of genes from the "known" (n = 180) senescence genes used in study
            hyp: hypergeometric p value for "known" gene overlap with hub genes
    '''
    
    
    def __init__(self, species = None, sig_type = 'hubs'):
        
        if species is None or (species != 'Human' and species != 'Mouse'):
            raise ValueError("Please pick Human or Mouse for species")

        if species == 'Mouse':
            #literature markers
            stream = resource_filename(__name__, 'data/Mouse_literature_markers.pickle')
            with open(stream, 'rb') as handle:
                self.literature_markers = pickle.load(handle) 
            #senGPT markers
            stream = resource_filename(__name__, 'data/senGPT.txt')
            with open(stream, 'r') as f:
                sengpt = [x.strip() for x in list(f)]
                self.senGPT = [x[0] + x[1:].lower() for x in sengpt]

    
        if species == 'Human':
            #literature markers
            stream = resource_filename(__name__, 'data/Human_literature_markers.pickle')
            with open(stream, 'rb') as handle:
                self.literature_markers = pickle.load(handle)
            #senGPT markers
            stream = resource_filename(__name__, 'data/senGPT.txt')
            with open(stream, 'r') as f:
                self.senGPT = [x.strip() for x in list(f)]
            

        #load Senpy Data
        if species == 'Mouse' and sig_type == 'hubs':
            stream = resource_filename(__name__, 'data/5_TMS_HUBS_DICTIONARY_FILTERED.pickle')
            with open(stream, 'rb') as handle:
                self.hubs = pickle.load(handle)
            stream = resource_filename(__name__, 'data/5_TMS_HUBS_METADATA_FILTERED.pickle')
            self.metadata = read_pickle(stream)
        if species == 'Mouse' and sig_type == 'cell_type':
            stream = resource_filename(__name__, 'data/5_TMS_SIGS_DICTIONARY_FILTERED.pickle')
            with open(stream, 'rb') as handle:
                self.hubs = pickle.load(handle)
            stream = resource_filename(__name__, 'data/5_TMS_SIGS_METADATA_FILTERED.pickle')
            self.metadata = read_pickle(stream)
        if species == 'Human' and sig_type == 'hubs':
            stream = resource_filename(__name__, 'data/6_HUMAN_HUBS_DICTIONARY_FILTERED.pickle')
            with open(stream, 'rb') as handle:
                self.hubs = pickle.load(handle)
            stream = resource_filename(__name__, 'data/6_HUMAN_HUBS_METADATA_FILTERED.pickle')
            self.metadata = read_pickle(stream)
        if species == 'Human' and sig_type == 'cell_type':
            stream = resource_filename(__name__, 'data/6_HUMAN_SIGS_DICTIONARY_FILT0255ERED.pickle')
            with open(stream, 'rb') as handle:
                self.hubs = pickle.load(handle)
            stream = resource_filename(__name__, 'data/6_HUMAN_SIGS_METADATA_FILTERED.pickle')
            self.metadata = read_pickle(stream)



    def get_genes(self, search_tuple):
        '''
        Gets a list of genes from the hub or signature

        search_tuple: tuple
            lookup tuple for your hub or signature of intrest.
            if load_hubs sig_type was 'hubs', then the format is (tissue, cell, hub_num) - (str, str, int)
            if load_hubs sig_type was 'cell_type', then the format is (tissue, cell) - (str, str)
        '''

        dict_format_hub = self.hubs[search_tuple]

        return [x[0] for x in dict_format_hub]


    
    def merge_hubs(self, meta_df, new_name = None,
                   overlap_threshold = 1,
                   calculate_thresh = False,
                   bg_N = 25000,
                   n_permutations = 1000,
                   p_thres = 0.05):
        '''
        Merge signatures together and keep only the genes that occur in overlap_threshold number of individual
        signatures. You can also cacluate the threshold based on random permutations with calculate_thresh
        set to True. The values in the resulting dictionary are the number of occurences

        meta_df: pandas dataframe from self.metadata
            Signatures to be merged. Each signature provided by each row in the metadata will be merged.
            Rows can be filtered, or passing the whole metadata will provide a universal signature.

        new_name: str
            What to name the new merged signature. It will be saved in self.hubs under this key name

        overlap_threshold: int
            How many of the individual signatures must a gene be present in to be kept in the merged
            signature. The default of 1 takes the union of all signatures. Ignored if calculate_thresh
            is set to True.

        calculate_thresh: bool
            Whether to calculate a threshold automatically using a permutation approach. This will find
            the number of signature occurences genes will have by random chance given the signatures
            provided. More useful when merging a higher number of signatures.

        bg_N: int
            Number of background genes to sample if using calculate_thresh.

        n_permutations: int
            Number of permutations used in calculate_thresh.

        p_thresh: float
            Value used to determine calculate_thres. The default of 0.05 (5%) means that the all genes
            retained occur more often than what is expected by at most a 5% chance. Actual chance may
            be smaller and is printed.
            
        '''
        
        if new_name is None:
            raise ValueError("Please provide a signature name via new_name")

        
        union_genes = []
        signatures = {}
        for x in range(len(meta_df)):
            tissue = meta_df.iloc[x].tissue
            cell = meta_df.iloc[x].cell
            if 'hub_num' in meta_df.columns:
                hub_num =  meta_df.iloc[x].hub_num
                union_genes += [x[0] for x in self.hubs[(tissue, cell, hub_num)]]
                signatures[(tissue, cell, hub_num)] = self.hubs[(tissue, cell, hub_num)]
            else:
                union_genes += [x[0] for x in self.hubs[(tissue, cell)]]
                signatures[(tissue, cell)] = self.hubs[(tissue, cell)]
         
            gene_counts = Counter(union_genes)

        #if calculate_thresh is False you return just a 
        if not calculate_thresh:
            gene_counts = {item: count for item, count in gene_counts.items() if count >= overlap_threshold}
            self.hubs[new_name] = list(gene_counts.items())


        if calculate_thresh:
            bg = range(bg_N)
            
            #map genes in sigs to arbitrary ids
            _ = list(set(union_genes))
            gene_ids = dict(zip(_, range(len(_))))

            all_genes_redundant = []
            for x in list(signatures):
                all_genes_redundant += [x[0] for x in signatures[x]]
            all_genes_redundant = [gene_ids[x] for x in all_genes_redundant]

            out = []
            for x in tqdm(range(n_permutations)):
                all_genes_redundant_temp = []
                for y in range(len(meta_df)):
                    s = meta_df.iloc[y].size
                    l = choices(bg, k = s)
                    all_genes_redundant_temp += l
                        
                all_genes_nr_temp = list(set(all_genes_redundant_temp))
                
            
                for gene in all_genes_nr_temp:
                    out += [all_genes_redundant.count(gene)]
                    
            out = array(out)

            for calculated_thresh in range(100):
                res = len(out[out >= calculated_thresh])/len(out)
                if res < p_thres:
                    break
            print(f'A gene will occur {calculated_thresh} times at {round(res*100,2)}% chance')
            print(f'Threfore {calculated_thresh} is the calculated_threshold')
            gene_counts = {item: count for item, count in gene_counts.items() if count >= calculated_thresh}
            self.hubs[new_name] = list(gene_counts.items())
        


    
    def search_hubs_by_genes(self, genes, bg_N = 25000):
        '''
        Pass a list of genes and return a dataframe of ranked hubs by overlap

        p_value is hypergoemetric survivial. p_adj is a BH correction.

        bg_N: int
            Number of genes in the background. Change to better match your context for
            a more accurate Hypergeometric P-value.
        
        '''
        
        out = []
        for hub in self.hubs:
            hub_genes = [x[0].lower() for x in self.hubs[hub]]
            
            size = len(hub_genes)
            
            overlap = [x for x in genes if x.lower() in hub_genes]

            sf = hypergeom_sf(len(overlap) - 1, bg_N, len(genes), size)
            
            out.append(list(hub) + [size, len(overlap), sf, overlap])
            
            df = DataFrame(out, columns = ['tissue', 'cell_type', 'hub_num' , 'size', 'num_hits', 'p_value', 'hits'])
            
            df = df.sort_values('p_value').reset_index(drop = True)
            df['p_adj'] = (df['p_value'] * len(df) / (df.index + 1)).clip(upper=1)

            df = df[['tissue', 'cell_type', 'hub_num' , 'size', 'num_hits', 'p_value', 'p_adj', 'hits']]

        return df
    