
import multiprocessing
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import fisher_exact
from joblib import Parallel, delayed


def rank_cluster_specific_variants(heteroplasmy_df, binary_cell_label, hetero_threshold=0.1, ncores=-1, 
	min_cell_number=6, flag_control_cell_number=5, keep_negative_enrichment=False, transition_only=True):
	"""Detect variants specific to a group of cells defined by binary_cell_label.
	Binarize the heteroplasmy df and run fisher's exact test for each variant.
	Note the function automatically takes the intersection of cells between heteroplasmy_df and binary_cell_label.
	:param heteroplasmy_df: Cell x variant matrix.
	:param binary_cell_label: Binary pd.Series with 1 denoting cells in group of interest and 0 comparison group
	:param hetero_threshold: Heteroplasmy threshold to use for binarization
	:param ncores: Number of cores to use for parallelized fisher's exact test across variants.
	:return: pd.DataFrame with odds ratio and p-val for all enriched variants (odds ratio > 1)
	"""
	if ncores == -1:
		ncores = multiprocessing.cpu_count()
	
	binary_cell_label = binary_cell_label[binary_cell_label.index.intersection(heteroplasmy_df.index)]
	print('{} cells in target group.'.format(binary_cell_label.sum()))
	print('{} cells in control group.'.format((~binary_cell_label).sum()))

	bin_hetero_df = (heteroplasmy_df > hetero_threshold).astype(int)
	bin_hetero_df = bin_hetero_df.loc[binary_cell_label.index, :]
	bin_hetero_df = bin_hetero_df.loc[:, bin_hetero_df.sum() > 0]
	
	# keep only variants present in less than 10% of all cells, otherwise odds ratio becomes unreliable
	frac_group = bin_hetero_df.sum() / bin_hetero_df.shape[0]
	frac_group = frac_group < 0.1
	bin_hetero_df = bin_hetero_df.loc[:, frac_group]

	# keep only variants present in at least min_cell_number of cells
	bin_hetero_df = bin_hetero_df.loc[:, bin_hetero_df.sum() >= min_cell_number]
	
	all_vars = bin_hetero_df.columns
	fisher_results = Parallel(n_jobs=ncores)(delayed(fisher_exact)(pd.crosstab(binary_cell_label, bin_hetero_df[var])) for var in all_vars)
	fisher_results = pd.DataFrame(fisher_results, index=bin_hetero_df.columns, columns=['odds_ratio', 'p_val'])
	fisher_results['target_count'] = bin_hetero_df.loc[binary_cell_label, fisher_results.index].sum()
	fisher_results['control_count'] = bin_hetero_df.loc[~binary_cell_label, fisher_results.index].sum()
	fisher_results = fisher_results.sort_values(['p_val'])

	# flag variants in too many control set
	fisher_results['flag'] = fisher_results['control_count'] > flag_control_cell_number
	
	# whether to keep variants with negative enrichment
	if not keep_negative_enrichment:
		fisher_results = fisher_results.loc[fisher_results['odds_ratio'] > 1, :]

	# keep only transition mutations, otherwise much more likely to be technical
	if transition_only:
		transition_mutations = ['A>G', 'G>A', 'T>C', 'C>T']
		var_use = pd.Index([x for x in fisher_results.index if any([y in x for y in transition_mutations])])
		fisher_results = fisher_results.loc[var_use, :]
	
	return fisher_results


def compute_affinity(heteroplasmy_df, dist='weighted_jaccard'):
	"""Compute affinity matrix using cell heteroplasmy df.
	:param heteroplasmy_df: Cell x variant matrix.
	:return: Cell-cell affinity_df.
	"""
	if dist == 'weighted_jaccard':
		dist = _weighted_jaccard

	affinity_df = squareform(pdist(heteroplasmy_df, dist))
	affinity_df = pd.DataFrame(affinity_df, index=heteroplasmy_df.index, columns=heteroplasmy_df.index)

	return affinity_df


def compute_cluster_clonality_index(affinity_df, cluster_assignment, scale_factor=10, round_digits=4):
	"""Quantify how clonal each cluster of cells is.
	Note the function automatically takes the intersection of cells between affinity_df and cluster_assignment.
	:param affinity_df: Cell-cell affinity_df.
	:param cluster_assignment: pd.Series encoding cluster assignment for cells.
	:return: A pd.Series encoding cluster clonality index for all clusters.
	"""
	cluster_clonality_index = dict()
	for i in sorted(cluster_assignment.unique()):
		use_index = cluster_assignment[cluster_assignment==i].index.intersection(affinity_df.index)
		cluster_clonality_index[i] = affinity_df.loc[use_index, use_index].sum().sum()  # sum of all edges within cluster
		cluster_clonality_index[i] /= len(use_index)**2  # normalize wrt number of nodes in the cluster
	
	cluster_clonality_index = (pd.Series(cluster_clonality_index) * scale_factor).round(round_digits)   # scale up
	cluster_clonality_index = cluster_clonality_index.sort_values(ascending=False)

	return cluster_clonality_index


def computer_intercluster_connectedness(affinity_df, cluster_assignment):
	return


def _weighted_jaccard(a, b):
	stacked = np.array([a, b])
	return (stacked.min(axis=0).sum() / stacked.max(axis=0).sum())