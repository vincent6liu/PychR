import pandas as pd
import numpy as np
import h5py
import scipy.sparse
import glob
from collections import defaultdict


def load_gs_mat_from_arrow(arrow_file_dir):
    """
    arrow_file_dir: directory containing .arrow files of all replicates/samples
    """
    if arrow_file_dir[-1] != '/':
        arrow_file_dir += '/'
    arrow_files = sorted(glob.glob(arrow_file_dir+'*.arrow'))
    # construct chr list in counting order, since it's alphanumeric by default
    chr_keys = ['chr'+str(i) for i in range(1, 23)] + ['chrX']
    combined_mat_list, cell_barcode_list = [], []
    
    for arrow_file_path in arrow_files:
        short_id = arrow_file_path.split('/')[-1].replace('.arrow', '')
        arrow_file = h5py.File(arrow_file_path, 'r')
        peak_mat = arrow_file['GeneScoreMatrix']
        
        cell_barcodes = [short_id+'#'+bc.decode() for bc in peak_mat['Info']['CellNames']]
        
        cell_num = len(cell_barcodes)
        chr_mat_list = []
        
        for ch in chr_keys:
            #print(ch)
            chr_mat_arrow = peak_mat[ch]
            x = chr_mat_arrow['x']
            i = chr_mat_arrow['i']
            dim = (len(chr_mat_arrow['rowSums'][0]), len(chr_mat_arrow['jLengths'][0]))
            idxptr = np.concatenate([[0], np.cumsum(chr_mat_arrow['jLengths'][0])])
            chr_peak_mat = scipy.sparse.csc_matrix((x[0, :], i[0, :]-1, idxptr), shape=dim).todense()
            
            if dim[1] != cell_num:
                # there a cells ommitted by ArchR due to all zero values for this chromosome
                jvals = chr_mat_arrow['jValues'][0]
                ommitted_idx = np.sort(list(set(range(1, cell_num+1)) - set(jvals)))
                
                for o_i in ommitted_idx:
                    if o_i != cell_num:
                        chr_peak_mat = np.insert(chr_peak_mat, o_i-1, np.zeros(dim[0]), axis=1)
                    else:
                        chr_peak_mat = np.append(chr_peak_mat, np.zeros((dim[0], 1)), axis=1)
            
            chr_mat_list.append(chr_peak_mat)
        
        combined_mat = np.vstack(chr_mat_list)
        combined_mat_list.append(combined_mat)
        cell_barcode_list.append(cell_barcodes)
    
    rep_merged_mat = np.hstack(combined_mat_list).T
    rep_merged_barcodes = [bc for rep_bc_list in cell_barcode_list for bc in rep_bc_list]
    chr_genes = [(x[0].decode(), x[4].decode()) for x in arrow_file['GeneScoreMatrix']['Info']['FeatureDF'][:]]
    chr_gene_dict = defaultdict(list)
    
    for chr_gene in chr_genes:
        chr_gene_dict[chr_gene[0]].append(chr_gene[1])
    
    genes = [gene for ch in chr_keys for gene in chr_gene_dict[ch]]
    
    return rep_merged_mat, rep_merged_barcodes, genes


def load_motif_mat_from_arrow(arrow_file_dir):
    """
    arrow_file_dir: directory containing .arrow files of all replicates/samples
    """
    if arrow_file_dir[-1] != '/':
        arrow_file_dir += '/'
    arrow_files = sorted(glob.glob(arrow_file_dir+'*.arrow'))
    deviation_mat_list, z_mat_list, cell_barcode_list = [], [], []
    
    for arrow_file_path in arrow_files:
        short_id = arrow_file_path.split('/')[-1].replace('.arrow', '')
        arrow_file = h5py.File(arrow_file_path, 'r')
        arrow_data = arrow_file['MotifMatrix']
        cell_barcodes = [short_id+'#'+bc.decode() for bc in arrow_data['Info']['CellNames']]
        cell_barcode_list.append(cell_barcodes)
        
        # get deviation matrix
        deviation_mat = arrow_data['deviations']
        x, i = deviation_mat['x'], deviation_mat['i']
        n_rows = deviation_mat['rowMeans'].shape[1]
        n_cols = int(deviation_mat['i'].shape[1]/n_rows)
        dim = (n_rows, n_cols)
        idxptr = np.concatenate([[0], np.cumsum(deviation_mat['jLengths'][0])])
        deviation_mat = scipy.sparse.csc_matrix((x[0, :], i[0, :]-1, idxptr), shape=dim).todense()
        deviation_mat_list.append(deviation_mat)
        
        # get zscore matrix
        z_features = [x[2].decode() for x in list(arrow_file['MotifMatrix']['Info']['FeatureDF']) if x[0].decode()=='z']
        z_mat = arrow_data['z']
        x, i = z_mat['x'], z_mat['i']
        n_rows = z_mat['rowMeans'].shape[1]
        n_cols = int(z_mat['i'].shape[1]/n_rows)
        dim = (n_rows, n_cols)
        idxptr = np.concatenate([[0], np.cumsum(z_mat['jLengths'][0])])
        z_mat = scipy.sparse.csc_matrix((x[0, :], i[0, :]-1, idxptr), shape=dim).todense()
        z_mat_list.append(z_mat)
        
    deviation_features = [x[2].decode() for x in list(arrow_data['Info']['FeatureDF']) if x[0].decode()=='deviations']
    z_features = [x[2].decode() for x in list(arrow_data['Info']['FeatureDF']) if x[0].decode()=='z']
    rep_merged_barcodes = [bc for rep_bc_list in cell_barcode_list for bc in rep_bc_list]
    
    rep_merged_deviation_mat = pd.DataFrame(np.hstack(deviation_mat_list).T, index=rep_merged_barcodes, columns=deviation_features)
    rep_merged_z_mat = pd.DataFrame(np.hstack(z_mat_list).T, index=rep_merged_barcodes, columns=deviation_features)
    
    return rep_merged_deviation_mat, rep_merged_z_mat


def load_qc_info_from_arrow(arrow_file_dir):
    """
    arrow_file_dir: directory containing .arrow files of all replicates/samples
    """
    if arrow_file_dir[-1] != '/':
        arrow_file_dir += '/'
    arrow_files = sorted(glob.glob(arrow_file_dir+'*.arrow'))
    # construct chr list in counting order, since it's alphanumeric by default
    meta_dict, cell_barcode_list = {}, []
    
    for arrow_file_path in arrow_files:
        short_id = arrow_file_path.split('/')[-1].replace('.arrow', '')
        arrow_file = h5py.File(arrow_file_path, 'r')
        meta_mat = arrow_file['Metadata']
        
        cell_barcodes = [short_id+'#'+bc.decode() for bc in meta_mat['CellNames']]
        
        doublet_enrichment = pd.Series(meta_mat['DoubletEnrichment'])
        doublet_score = pd.Series(meta_mat['DoubletScore'])
        tss_enrichment = pd.Series(meta_mat['TSSEnrichment'])
        n_frags = pd.Series(meta_mat['nFrags'])
        
        meta_df = pd.concat([doublet_enrichment, doublet_score, tss_enrichment, n_frags], axis=1)
        meta_df.index = cell_barcodes
        meta_df.columns = ['DoubletEnrichment', 'DoubletScore', 'TSSEnrichment', 'nFrags']
        
        meta_dict[short_id] = meta_df
    
    return meta_dict

