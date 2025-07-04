import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns
import numpy as np
import glob

def parseuropa(tsv):
    df = pd.read_table(tsv)
    df['peak'] = df['peak_chr'] + ':' + df['peak_start'].astype(str) + '-' + df['peak_end'].astype(str)
    df = df[
        ['peak', 'peak_chr', 'peak_start', 'peak_end', 'relative_location', 'distance', 'gene_id', 'gene_name', 'feature', 'transcript_biotype']
    ]
    # Parse uropa df to combine feature and rel_loc
    annFeature = []
    for i,r in df.iterrows():
        if r['relative_location'] in [
            'OverlapStart',
            'OverlapEnd',
            'PeakInsideFeature',
            'FeatureInsidePeak'
        ]:
            if r['relative_location'] == 'OverlapStart' and r['feature'] == 'transcript':
                annFeature.append('promoter')
            elif r['relative_location'] == 'OverlapStart' and r['feature'] == 'five_prime_utr':
                annFeature.append('promoter')
            else:
                annFeature.append(r['feature'])
        else:
            annFeature.append('intergenic')
    df['annFeature'] = annFeature
    return(df)

FLASHgroups = {
    'C8_MSL1': ['Dme_Clone8_MSL1_rep1', 'Dme_Clone8_MSL1_rep2', 'Dme_Clone8_MSL1_rep3'],
    'S2_MSL1': ['Dme_S2_MSL1_rep1', 'Dme_S2_MSL1_rep2', 'Dme_S2_MSL1_rep3'],
    'MSL2': ['Dme_MSL2_male_rep1', 'Dme_MSL2_male_rep2'],
    'MLE': ['Dme_MLE_male_rep1' , 'Dme_MLE_male_rep2']
}

bdir = snakemake.params.annodir
tdfcdic = {}
tdfdic = {}
for _g in FLASHgroups:
    dfs = []
    genes = []
    for _tdfname in FLASHgroups[_g]:
        _tdf = parseuropa(os.path.join(bdir, _tdfname + '_finalhits.txt'))
        _tdf['dataset'] = _tdfname
        dfs.append(_tdf)
        genes.append(set(_tdf[_tdf['annFeature'] != 'intergenic']['gene_name'].unique().tolist()))
    countdic = {}
    for _l in genes:
        for gene in _l:
            if gene not in countdic:
                countdic[gene] = 1
            else:
                countdic[gene] += 1
    gois = []
    for gene in countdic:
        if countdic[gene] == len(genes):
            gois.append(gene)
    fdf = pd.concat(dfs)
    ipname = '_'.join(_tdfname.split('_')[:-1])
    tdfcdic[ipname] = fdf[fdf['gene_name'].isin(gois)]['transcript_biotype'].value_counts()
    tdfdic[ipname] = fdf[fdf['gene_name'].isin(gois)]

# Get relative frequencies for annotations
FLASH_anno = pd.DataFrame(tdfcdic).fillna(0).transform(lambda x: x / x.sum()).T.reset_index().melt(id_vars='index')
FLASH_anno.columns = ['Sample', 'transcript_biotype', 'frequency']
# Get some simple metrics
mdic = {}
for ipname in tdfdic:
    mdic[ipname] = {}
    mdic[ipname]['CL_per_gene'] = tdfdic[ipname]['gene_name'].value_counts().mean()
    mdic[ipname]['CL'] = len(tdfdic[ipname])
    mdic[ipname]['genes'] = len(tdfdic[ipname]['gene_name'].unique())


fig = plt.figure(figsize=(8, 4), layout="constrained")
spec = fig.add_gridspec(2, 3)

ax0 = fig.add_subplot(spec[0, 0])
ax1 = fig.add_subplot(spec[0, 1])
ax2 = fig.add_subplot(spec[0, 2])
ax3 = fig.add_subplot(spec[1, :])

sns.barplot(
    data=pd.DataFrame(mdic).T.reset_index(),
    hue='index',
    y='CL_per_gene',
    ax=ax0,
    legend=False
)
ax0.set_title('Cross-links per gene')
ax0.set_ylabel('')

sns.barplot(
    data=pd.DataFrame(mdic).T.reset_index(),
    hue='index',
    y='CL',
    ax=ax1,
    legend=False
)
ax1.set_title('Total cross-links')
ax1.set_ylabel('')

sns.barplot(
    data=pd.DataFrame(mdic).T.reset_index(),
    hue='index',
    y='genes',
    ax=ax2
)
ax2.set_title('Total number of genes')
ax2.set_ylabel('')
ax2.legend(bbox_to_anchor=(1, 0.8))

sns.barplot(
    data=FLASH_anno,
    x='Sample',
    hue='transcript_biotype',
    y='frequency',
    ax=ax3,
    hue_order=['snoRNA', 'protein_coding', 'pseudogene', 'ncRNA', 'snRNA', 'pre_miRNA', 'tRNA']
)
ax3.legend(bbox_to_anchor=(1, 1))
ax3.set_xlabel('')
fig.savefig(snakemake.output.opng, dpi=300)