import pickle
import pandas as pd
import pkg_resources

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
            n_sen: number of genes from the "known" (n = 180)senescence genes used in study
            hyp: hypergeometric p value for "knwon" gene overlap with hub genes
    '''
    
    
    def __init__(self, species = None, sig_type = 'hubs'):
        
        if species is None or (species != 'Human' and species != 'Mouse'):
            raise ValueError("Please pick Human or Mouse for species")
        
        
        if species == 'Mouse' and sig_type == 'hubs':
            stream = pkg_resources.resource_filename(__name__, 'data/5_TMS_HUBS_DICTIONARY_FILTERED.pickle')
            with open(stream, 'rb') as handle:
                self.hubs = pickle.load(handle)
            stream = pkg_resources.resource_filename(__name__, 'data/5_TMS_HUBS_METADATA_FILTERED.pickle')
            self.metadata = pd.read_pickle(stream)
            
        if species == 'Mouse' and sig_type == 'cell_type':
            stream = pkg_resources.resource_filename(__name__, 'data/5_TMS_SIGS_DICTIONARY_FILTERED.pickle')
            with open(stream, 'rb') as handle:
                self.hubs = pickle.load(handle)
            stream = pkg_resources.resource_filename(__name__, 'data/5_TMS_SIGS_METADATA_FILTERED.pickle')
            self.metadata = pd.read_pickle(stream)
            
        if species == 'Human' and sig_type == 'hubs':
            stream = pkg_resources.resource_filename(__name__, 'data/6_HUMAN_HUBS_DICTIONARY_FILTERED.pickle')
            with open(stream, 'rb') as handle:
                self.hubs = pickle.load(handle)
            stream = pkg_resources.resource_filename(__name__, 'data/6_HUMAN_HUBS_METADATA_FILTERED.pickle')
            self.metadata = pd.read_pickle(stream)
            
        if species == 'Human' and sig_type == 'cell_type':
            stream = pkg_resources.resource_filename(__name__, 'data/6_HUMAN_SIGS_DICTIONARY_FILTERED.pickle')
            with open(stream, 'rb') as handle:
                self.hubs = pickle.load(handle)
            stream = pkg_resources.resource_filename(__name__, 'data/6_HUMAN_SIGS_METADATA_FILTERED.pickle')
            self.metadata = pd.read_pickle(stream)
        
    def search_hubs_by_genes(self, genes):
        '''
        Pass a list of genes and return a dataframe of ranked hubs by overlap
        
        '''
        
        out = []
        for hub in self.hubs:
            hub_genes = [x[0].lower() for x in self.hubs[hub]]
            
            size = len(hub_genes)
            
            overlap = [x for x in genes if x.lower() in hub_genes]
            
            out.append(list(hub) + [size, len(overlap), overlap])
            
            df = pd.DataFrame(out, columns = ['tissue', 'cell_type', 'hub_num' , 'size', 'num_hits', 'hits'])
            
            df = df.sort_values(['num_hits', 'size'], ascending = [False, True]).reset_index(drop = True)
            
        return df
    