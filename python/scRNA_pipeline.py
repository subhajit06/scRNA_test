import numpy as np
import pandas as pd
import scanpy as sc
sc.settings.verbosity = 3  
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80)

import numpy as np
import pandas as pd
import scanpy as sc
sc.pp.normalize_per_cell(adata, counts_per_cell_after = 1e4 )
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean = 0.0125 , max_mean = 3 , min_disp = 0.5)
sc.pp.regress_out(adata, ['n_counts']
sc.pp.scale(adata, max_value = 10 ) 
sc.tl.pca(adata, svd_solver = 'arpack' )
sc.pp.neighbors(adata, n_neighbors = 30 , n_pcs = 50 ) 
sc.tl.tsne(adata)
sc.tl.louvain(adata)
sc.tl.rank_genes_groups(adata,  'louvain' , method = 't-test' )
