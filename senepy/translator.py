import numpy as np
import pandas as pd
import pkg_resources



class translator(object):
    '''
    
    Mapper between ids in your dataset and hub genes. This is an automatic attempt that accounts
    for differences in capitalization and some aliases. You can manually add to the mapper through
    the self.add_map() function.
    
    hub: hub dictionary of multiple hubs or a single hub list
    
    data: addata object
        The hub gene names will be mapped to the gene names in these var_names
        
    '''
    def __init__(self, hub = None, data = None):

        
        if hub is not None:
            #break hub in to genes
            self.hub_genes = []
            
            if isinstance(hub, dict):
                for hu in hub:
                    self.hub_genes += [x[0] for x in hub[hu]]
            elif isinstance(hub, list):
                self.hub_genes = [x[0] for x in hub]
            else:
                raise ValueError("hub needs to be a hub dictionary or a single hub list")
                
                
        self.initial_genes_not_in_data = [x for x in self.hub_genes if x not in data.var_names]
        
        print(f'{len(self.initial_genes_not_in_data)} of {len(self.hub_genes)} genes not initially present')
        
        
        self.var_names_lower  = data.var_names.map(lambda x: x.lower())
        
        
        
        stream = pkg_resources.resource_stream(__name__, 'data/mart_export.txt')
        self.aliases = pd.read_csv(stream)
        self.aliases['Gene Synonym'] = self.aliases['Gene Synonym'].str.lower()
        self.aliases['Gene name'] = self.aliases['Gene name'].str.lower()
        self.aliases['UniProtKB Gene Name symbol'] = self.aliases['UniProtKB Gene Name symbol'].str.lower()

        self.aliases = self.aliases[(self.aliases['Gene Synonym'].isin(self.var_names_lower))\
            | (self.aliases['Gene name'].isin(self.var_names_lower))]

        
        
        self.aliases_d = dict(zip(self.aliases['Gene Synonym'], self.aliases['Gene name']))
        self.aliases_rev_d = dict(zip(self.aliases['Gene name'], self.aliases['Gene Synonym']))
        
        
        
        def translate_gene(var_names, var_names_lower, gene, aliases_d, aliases_rev_d):
    
            gene_lower = gene.lower()

            if gene_lower in var_names_lower:
                i = np.where(var_names_lower == gene.lower())[0][0]
                return var_names[i]

            try:
                alias =  aliases_d[gene_lower]
                if alias in var_names_lower:
                    i = np.where(var_names_lower == alias)[0][0]
                    return var_names[i]
            except:
                pass


            try:
                alias =  aliases_rev_d[gene_lower]
                if alias in var_names_lower:
                    i = np.where(var_names_lower == alias)[0][0]
                    return var_names[i]
            except:
                pass


            return 'not there'
        
        
        
        self.translated_genes = []
        self.mapper = {}
        self.untranslated_genes = []
        for gene in self.initial_genes_not_in_data:
            
            mapped = translate_gene(data.var_names, self.var_names_lower, gene, self.aliases_d, self.aliases_rev_d)
            
            self.translated_genes.append(mapped)
            
            if mapped != 'not there':
                self.mapper[gene] = mapped
            else:
                self.untranslated_genes.append(gene)
            
        self.translated_genes = [x for x in self.translated_genes if x != 'not there']
        
            
        print(f'{len(self.translated_genes)} of {len(self.initial_genes_not_in_data)} translated')
        print(f'{len(self.initial_genes_not_in_data)- len(self.translated_genes)} still not present')
        
        
    def add_map(self, d):
        '''
        manually pass a dictionary to appened to self.mapper
        '''
        self.mapper = self.mapper | d
            









