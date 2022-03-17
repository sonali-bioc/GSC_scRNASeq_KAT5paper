#------------------------------------------------------
# RNA velocity analysis - scVelo
#-----------------------------------------------------

# for KAT5 Project
# Samantha O'Connor
# Updated 03/17/22
# Example: GSC827_inVitro_combined

#------------------------------------------------------
# Setup
#---------------------------------------------------

# Files need before begin: velocyto loom file, loom file from Seurat completed analysis

docker run -it -v '/home/soconnor/:/files' cplaisier/scrna_seq_velocity_6_7_2021
cd /files
pip3 install -U umap-learn==0.3.10
python3

# Imports
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import loompy


#------------------------------------------------------
# Load files
#---------------------------------------------------
dir1 = 'CytmRNA-CD8KO_onlyhg38'
dir2 = 'CytmRNA-KAT5KO_onlyhg38'
data = 'GSC827_inVitro_combined'

adata_Seurat = scv.read_loom(dir2+'/'+data+'_data_POST.loom')
umap_coord = pd.read_csv(dir2+'/'+data+'_UMAPCoordinates.csv', header = 0, index_col = 0) # if want to use pre-defined coordinates
adata_Seurat.obsm['umap_cell_embeddings'] = np.array(umap_coord.loc[adata_Seurat.obs_names])
​
#---- Run datasets (CD8KO, KAT5KO) independently: -----
# Load data from velocyto
fold1 = 'velocyto'
#adata_vc_CD8KO = scv.read_loom(dir1+'/'+fold1+'/'+dir1+'.loom')
adata_vc_KAT5KO = scv.read_loom(dir2+'/'+fold1+'/'+dir2+'.loom')
adata_vc = adata_vc_KAT5KO

#---- Run datasets (CD8KO, KAT5KO) together: -----
# Merge two velocyto loom files
os.chdir(dir2)
files = ["CD8KO.loom","KAT5KO.loom"] # copied/renamed and put them in same folder
loompy.combine(files, "merged.loom", key="Accession")
# OR
# on the command line in folder do: cp file1.loom merged.loom
#ds = loompy.connect("merged.loom")
#for fn in files[1:]:
#    ds.add_loom(fn, batch_size=1000)

adata_vc = scv.read_loom('merged.loom')

#------------------------------------------------------
# Merge data
#---------------------------------------------------
tag = "GSC827_inVitro_combined"

adata = scv.utils.merge(adata_vc, adata_Seurat)
adata.obs["seurat_new_clusters"]=adata.obs["seurat_clusters"]-1
adata.obs.seurat_new_clusters = adata.obs.seurat_new_clusters.astype('category')

# Colorize by UMAP colors
os.chdir('..')
dir3 = 'scRNA_seq_Paddison/redo_analysis_with_ASU_filters'
ident_colors = pd.read_csv(dir3+'/'+tag+'_seurat1_seurat_umap_colors.csv')
adata.uns['ClusterName_colors']=list(ident_colors["ucols"])


#------------------------------------------------------
# Preprocess the data
#---------------------------------------------------

scv.pp.filter_and_normalize(adata, min_shared_counts=30, n_top_genes=3000)
scv.pp.moments(adata)
​
#------------------------------------------------------
# Estimate RNA velocity
#---------------------------------------------------

os.chdir(tag)
newpath = 'driverGenes'
if not os.path.exists(newpath):
    os.makedirs(newpath)

scv.tl.recover_dynamics(adata, n_jobs=8) # run by itself

modes = ['dynamical', 'stochastic']
for mode1 in modes:
    scv.tl.velocity(adata, mode=mode1)
    scv.tl.velocity_graph(adata, basis='umap_cell_embeddings')
    # Plot velocity with streamlines
    scv.pl.velocity_embedding_stream(adata, basis='umap_cell_embeddings', color=['seurat_new_clusters'], palette = adata.uns['ClusterName_colors'], save = tag+'_velocity_embedding_'+mode1+'.png', dpi=300)
    # Plot velocity with arrows
    scv.pl.velocity_embedding(adata, basis='umap_cell_embeddings', color=['seurat_new_clusters'], palette = adata.uns['ClusterName_colors'], arrow_length=5, arrow_size=3, dpi=300, save = tag+'_velocity_embedding_stream_arrows_'+mode1+'.png')
    # Identify velocity driver genes
    scv.tl.rank_velocity_genes(adata, groupby='seurat_new_clusters', min_corr=.3)
    df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
    df.to_csv(newpath+'/'+tag+'_rank_velocity_genes_'+mode1+'.csv')
    for clust in list(range(0, len(adata.obs['seurat_new_clusters'].value_counts()))):
        n = 10 # plot top 10 velocity genes
        scv.pl.velocity(adata, list(df[str(clust)][0:n]), ncols=2, add_outline=True, color=['seurat_new_clusters'], palette = adata.uns['ClusterName_colors'], basis='umap_cell_embeddings', save = tag+ '_velocity_genes_cluster_'+str(clust)+'_top_'+str(n)+'_'+mode1+'.pdf')




​
