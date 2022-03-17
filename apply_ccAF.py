
#------------------------------------------------------
# ccAF application
#-----------------------------------------------------

# for KAT5 Project
# Samantha O'Connor
# Updated 03/17/22
# ccAF description and build information: O'Connor et al., 2021
# Example: GSC827_inVitro_combined


#------------------------------------------------------
# Setup
#----------------------------------------------------

# File needed before begin: pre-normalized scRNA-seq data with genes as ensembl IDs as LOOM file.

docker pull cplaisier/scrna_seq_velocity # if haven't already pulled
# Run the docker container using the following command (replace with home directory / parent folder of folder where scRNA-seq loom file is):
docker run -it -v '/home/soconnor:/files' cplaisier/scrna_seq_velocity
cd /files
python3

# Imports
import pandas as pd
import scanpy as sc
import ccAF


#------------------------------------------------------
# Read in file and apply ccAF
#---------------------------------------------------

dir1 = 'scRNA_seq_Paddison/redo_analysis_with_ASU_filters'
data1 = 'GSC827_inVitro_combined'

# Load in data
data = sc.read_loom(dir1+'/'+data1+'_data_check.loom')
data

# Apply ccAF
data.obs['ccAF'] = ccAF.ccAF.predict_labels(data)
data.obs['ccAF'].value_counts()
pd.DataFrame(data.obs['ccAF']).to_csv(dir1+'/'+data1+'_ccAF_calls_check.csv')
