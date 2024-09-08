import pandas as pd
import os
import numpy as np
import multiprocessing
import sklearn
from itertools import repeat
import pkg_resources
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i','--input_folder', action='store', dest='input_folder', help='the path contains the comparison stats files')
parser.add_argument('-i','--mode', action='store', dest='input_folder', help='the path contains the comparison stats files')
parser.add_argument('-a','--comp1_file', action='store', dest='comp1_file', help='first comparison file')
parser.add_argument('-b','--comp2_file', action='store', dest='comp2_file', help='second comparison file')
parser.add_argument('-o','--output_file', action='store', dest='output_file', help='output file, each row assigned with a c-score, a p-value, and a sense field. If the mode is gene, then also with a pc filed indicating whether a gene (row) is protein-coding.')
parser.add_argument("-m", "--mode", choices=['gene,pathway'], dest='mode', default ='gene')

paras = parser.parse_args()

comp1_file = paras.comp1_file
comp2_file = paras.comp2_file
output_file = paras.output_file
mode = paras.mode

# make a summary list of coding genes
# in the hope of increasing the significance of the pathways

DATA_PATH = pkg_resources.resource_filename('cscore', 'data/')
GTF_FILE = pkg_resources.resource_filename('cscore', 'data/gencode_v42.gtf')

gencode = pd.read_table(GTF_FILE, comment="#",
                        sep="\t", names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])

gencode_genes = gencode[(gencode.feature == "transcript")][[
    'seqname', 'start', 'end', 'attribute']].copy().reset_index().drop('index', axis=1)


def gene_info(x):
    # Extract gene names, gene_type, gene_status and level
    g_name = list(filter(lambda x: 'gene_name' in x,  x.split(";")))[
        0].split(" ")[2].strip('\"')
    g_type = list(filter(lambda x: 'gene_type' in x,  x.split(";")))[
        0].split(" ")[2]
    return (g_name, g_type)


gencode_genes["gene_name"], gencode_genes["gene_type"] = zip(
    *gencode_genes.attribute.apply(lambda x: gene_info(x)))
pc_genes = gencode_genes.query(
    "gene_type=='\"protein_coding\"'")
pc_gene_set = set(pc_genes.gene_name)


def weight(fdr):
    return np.where(fdr < 0.05, 1, np.log10(fdr)/np.log10(0.05))
def ratio(fc_comp1, fc_comp2):
    df = pd.DataFrame({'fc_comp1': fc_comp1, 'fc_comp2': fc_comp2})
    return np.where(fc_comp1*fc_comp2 > 0, np.maximum(np.abs(fc_comp1), np.abs(fc_comp2))/(np.abs(fc_comp1-fc_comp2)+1), -np.abs(fc_comp1-fc_comp2)/(np.maximum(np.abs(fc_comp1), np.abs(fc_comp2))+1))
def score(comp1_np, comp2_np):
    fc_comp1 = comp1_np[:,0]
    fc_comp2 = comp2_np[:,0]
    fdr_comp1_weight = weight(comp1_np[:,1])
    fdr_comp2_weight = weight(comp2_np[:,1])
    magnitude = np.abs(fc_comp1*fdr_comp1_weight) + np.abs(fdr_comp2_weight*fc_comp2)
    return np.array(magnitude * ratio(fc_comp1, fc_comp2))

def shuffle_calc(step, comp1_np, comp2_np):
    comp1_shuffle = sklearn.utils.shuffle(comp1_np, random_state=step)
    comp2_shuffle = sklearn.utils.shuffle(comp2_np, random_state=step+40000)
    permutation_score = score(comp1_shuffle, comp2_shuffle)
    return permutation_score

if __name__ == '__main__':
    comp1 = pd.read_csv(comp1_file, sep='\t')
    comp2 = pd.read_csv(comp2_file, sep='\t')
    rows_keep = np.intersect1d(comp1['Unnamed: 0'], comp2['Unnamed: 0'])
## same genes and same order
    comp1 = comp1[comp1['Unnamed: 0'].isin(rows_keep)]
    comp2 = comp2[comp2['Unnamed: 0'].isin(rows_keep)]
## order dataframes by Unnamed: 0
    comp1 = comp1.sort_values('Unnamed: 0').reset_index(drop=True)
    comp2 = comp2.sort_values('Unnamed: 0').reset_index(drop=True)
    comp1['avg_log2FC'].astype(float)
    comp2['avg_log2FC'].astype(float)
    comp1['p_val_adj'].astype(float)
    comp2['p_val_adj'].astype(float)
    comp1_np = comp1[['avg_log2FC', 'p_val_adj']].to_numpy(dtype=float)
    comp2_np = comp2[['avg_log2FC', 'p_val_adj']].to_numpy(dtype=float)
    scores = np.array(score(comp1_np, comp2_np))
        ## filter with score not equal to 0
        ## filter comp1 and comp2 dataframe rows
        # filter_cond = (scores!=0) & ((comp1['p_val_adj']<0.05) | (comp2['p_val_adj']<0.05))
    filter_cond = scores!=0
        # print(sum(filter_cond))
    comp1 = comp1.loc[filter_cond]
    comp2 = comp2.loc[filter_cond]
    comp1_np = comp1_np[filter_cond]
    comp2_np = comp2_np[filter_cond]
    scores_notzero = scores[filter_cond]
        # scores_notzero = scores
## calculate permutatoin scores
    comp1_shuffles = []
    comp2_shuffles = []
## TODO: can also rewrite to yield
## if len(genes) < 200, k = len(genes)**2, else k = 40000
    if comp1_np.shape[0] < 200:
        k = comp1_np.shape[0]**2
    else:
        k = 40000
    pocomp2 = multiprocessing.Pocomp2(64)
    permutation_scores = np.vstack(pocomp2.starmap(shuffle_calc, zip(range(0, int(k)), repeat(comp1_np), repeat(comp2_np))))
    scores_all = np.vstack((scores_notzero, permutation_scores))
## the p-value is the proportion of the permutation scores that are more extreme than the original score
## these fcomp2lowing calculation should be row-wise
## for the scores array, if the first ccomp2umn is positive, then the p-value is the proportion of the permutation scores that are greater than the original score
## if the first ccomp2umn is negative, then the p-value is the proportion of the permutation scores that are less than the original score
    ps = np.ones(scores_all.shape[1], dtype = float)
    sense = []
    for i in range(scores_all.shape[1]):
        ps_tmp = float(np.sum(scores_all[1:, i] > scores_all[0, i]))/float((scores_all.shape[0]-1))
        if ps_tmp < 0.5:
            sense.append('high')
            ps[i] = ps_tmp
        else:
            sense.append('low')
            ps[i] = 1-ps_tmp
    score_p = np.vstack((scores_all[0, :], ps))
## merge comp1 and comp2 by Unnamed: 0
    df_comp1_comp2 = pd.merge(comp1, comp2, on='Unnamed: 0', suffixes=('_comp1', '_comp2'))
    df_comp1_comp2 = df_comp1_comp2.sort_values('Unnamed: 0').reset_index(drop=True)
    df_comp1_comp2['score'] = score_p[0, :]
    df_comp1_comp2['p'] = score_p[1, :]
    df_comp1_comp2['convergence'] = sense
## get annotation whether the gene is coding or non-coding
## df where Unnamed: 0 is in pc_gene_set, the coding is True, else False
    if mode == 'gene':
        df_comp1_comp2['coding'] = df_comp1_comp2['Unnamed: 0'].isin(pc_gene_set)

    df_comp1_comp2.to_csv(output_file, index=False, float_format="%.64f", sep='\t')