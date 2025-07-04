import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os
import seaborn as sns
import numpy as np
import glob

cdic = {
    "three_prime_utr":"#ff7f0e",
    "five_prime_utr":"#9467bd",
    "first_exon":"#2ca02c",
    "intergenic":"#7f7f7f",
    "promoter":"#1f77b4",
    "first_intron":"#d62728",
    "other_exon": "#8c564b",
    "other_intron": "#e377c2"
}


def parsegtf(gtf):
    '''
    create genedic
    genedic[geneid] = {
        txn = (start, end)
        exons = [(start, end), (start, end), ...]
    }
    '''
    genedic = {}
    with open(gtf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            geneid = None
            l = line.strip().split('\t')
            for i in l[8].split(';'):
                if 'gene_id' in i:
                    geneid = i.split(' ')[-1].replace('"', '').strip()
            if geneid:
                if geneid not in genedic:
                    genedic[geneid] = {
                        'txn': None,
                        'exons': []
                    }
                if l[2] == 'transcript':
                    start = int(l[3])
                    end = int(l[4])
                    if not genedic[geneid]['txn']:
                        genedic[geneid]['txn'] = (start, end)
                    else:
                        if start < genedic[geneid]['txn'][0]:
                            genedic[geneid]['txn'] = (start, genedic[geneid]['txn'][1])
                        if end > genedic[geneid]['txn'][1]:
                            genedic[geneid]['txn'] = (genedic[geneid]['txn'][0], end)
                elif l[2] == 'exon':
                    start = int(l[3])
                    end = int(l[4])
                    genedic[geneid]['exons'].append((start, end))
    for i in genedic:
        genedic[i]['exons'] = sorted(genedic[i]['exons'], key=lambda x: x[0])
        # define first intron
        if len(genedic[i]['exons']) > 1:
            genedic[i]['firstintron'] = (genedic[i]['exons'][0][1], genedic[i]['exons'][1][0])
    return genedic

def parseuropa(tsv, genedic):
    df = pd.read_table(tsv)
    df['peak'] = df['peak_chr'] + ':' + df['peak_start'].astype(str) + '-' + df['peak_end'].astype(str)
    df = df[
        ['peak', 'peak_chr', 'peak_start', 'peak_end', 'relative_location', 'gene_id', 'gene_name', 'feature', 'transcript_biotype']
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
            annofeat = None
            if r['relative_location'] == 'OverlapStart' and r['feature'] == 'transcript':
                annFeature.append('promoter')
            elif r['relative_location'] == 'OverlapStart' and r['feature'] == 'five_prime_utr':
                annFeature.append('promoter')
            else:
                # Differentiate between exons and introns, but only if feature is transcript or exon.
                if r['feature'] in ['transcript', 'exon']:
                    geneid = r['gene_id']
                    assert geneid in genedic, f"Gene {geneid} not found in genedic.."
                    start = r['peak_start']
                    end = r['peak_end']
                    if r['feature'] == 'transcript':
                        # Since feature is transcript, it should not be in an exon.
                        for exon in genedic[geneid]['exons']:
                            assert start <= exon[1] or start >= exon[1], f"Feature transcript, but found in exon: {r} <-> {exon}.."
                        assert 'firstintron' in genedic[geneid], f"No first intron annotated, but feature is transcript.. {r} <-> {genedic[geneid]}"
                        if start >= genedic[geneid]['firstintron'][0] and start <= genedic[geneid]['firstintron'][1]:
                            annofeat = 'first_intron'
                        elif end >= genedic[geneid]['firstintron'][0] and end <= genedic[geneid]['firstintron'][1]:
                            annofeat = 'first_intron'
                        else:
                            annofeat = 'other_intron'
                    elif r['feature'] == 'exon':
                        if start >= genedic[geneid]['exons'][0][0] and start <= genedic[geneid]['exons'][0][1]:
                            annofeat = 'first_exon'
                        elif end >= genedic[geneid]['exons'][0][0] and end <= genedic[geneid]['exons'][0][1]:
                            annofeat = 'first_exon'
                        else:
                            annofeat = 'other_exon'
                    annFeature.append(annofeat)
                else:
                    annFeature.append(r['feature'])
        else:
            annFeature.append('intergenic')
    df['annFeature'] = annFeature
    return(df)

def plotter(MSL1, MSL1_RNA, padj_cutoff, genes, ofile):
    genedic = parsegtf(genes)

    MSL1 = parseuropa(MSL1, genedic)
    MSL1_RNA = pd.read_table(MSL1_RNA)

    # Define lists
    MSL1_RNA_downgenes = list(MSL1_RNA[(MSL1_RNA['padj'] < padj_cutoff) & (MSL1_RNA['log2FoldChange'] < 0)].index)
    downgenes_withpeaks = MSL1[MSL1['gene_id'].isin(MSL1_RNA_downgenes)]['gene_id'].unique()
    downgenes_withoutpeaks = [i for i in MSL1_RNA_downgenes if i not in downgenes_withpeaks]
    assert(len(downgenes_withpeaks) + len(downgenes_withoutpeaks) == len(MSL1_RNA_downgenes))

    downgenes_withpeaks_X = MSL1[ (MSL1['gene_id'].isin(downgenes_withpeaks)) & (MSL1['peak_chr'] == 'X') ]['gene_id'].unique()
    downgenes_withpeaks_nonX = MSL1[ (MSL1['gene_id'].isin(downgenes_withpeaks)) & (MSL1['peak_chr'] != 'X') ]['gene_id'].unique()
    assert(len(downgenes_withpeaks_X) + len(downgenes_withpeaks_nonX) == len(downgenes_withpeaks))

    downgenes_withoutpeaks_X = list(MSL1_RNA.loc[downgenes_withoutpeaks][MSL1_RNA.loc[downgenes_withoutpeaks]['Chr'] == 'X'].index)
    downgenes_withoutpeaks_nonX = list(MSL1_RNA.loc[downgenes_withoutpeaks][MSL1_RNA.loc[downgenes_withoutpeaks]['Chr'] != 'X'].index)
    assert(len(downgenes_withoutpeaks_X) + len(downgenes_withoutpeaks_nonX) == len(downgenes_withoutpeaks))

    # Define lists - upregulated genes
    MSL1_RNA_upgenes = list(MSL1_RNA[(MSL1_RNA['padj'] < padj_cutoff) & (MSL1_RNA['log2FoldChange'] > 0)].index)
    upgenes_withpeaks = MSL1[MSL1['gene_id'].isin(MSL1_RNA_upgenes)]['gene_id'].unique()
    upgenes_withoutpeaks = [i for i in MSL1_RNA_upgenes if i not in upgenes_withpeaks]
    assert(len(upgenes_withpeaks) + len(upgenes_withoutpeaks) == len(MSL1_RNA_upgenes))

    upgenes_withpeaks_X = MSL1[ (MSL1['gene_id'].isin(upgenes_withpeaks)) & (MSL1['peak_chr'] == 'X') ]['gene_id'].unique()
    upgenes_withpeaks_nonX = MSL1[ (MSL1['gene_id'].isin(upgenes_withpeaks)) & (MSL1['peak_chr'] != 'X') ]['gene_id'].unique()
    assert(len(upgenes_withpeaks_X) + len(upgenes_withpeaks_nonX) == len(upgenes_withpeaks))

    upgenes_withoutpeaks_X = list(MSL1_RNA.loc[upgenes_withoutpeaks][MSL1_RNA.loc[upgenes_withoutpeaks]['Chr'] == 'X'].index)
    upgenes_withoutpeaks_nonX = list(MSL1_RNA.loc[upgenes_withoutpeaks][MSL1_RNA.loc[upgenes_withoutpeaks]['Chr'] != 'X'].index)
    assert(len(upgenes_withoutpeaks_X) + len(upgenes_withoutpeaks_nonX) == len(upgenes_withoutpeaks))

    # Plot figure

    counts = np.array([
        [len(downgenes_withpeaks_X), len(downgenes_withpeaks_nonX)],
        [len(downgenes_withoutpeaks_nonX), len(downgenes_withoutpeaks_X)]
    ])

    counts_up = np.array([
        [len(upgenes_withpeaks_X), len(upgenes_withpeaks_nonX)],
        [len(upgenes_withoutpeaks_nonX), len(upgenes_withoutpeaks_X)]
    ])

    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(12,8))

    ######################################## Pie plots ########################################

    ax[0, 0].pie(
        counts.flatten(), radius=1-0.4,
        wedgeprops=dict(width=0.4, edgecolor='w'),
        colors=['#1f77b490', '#1f77b450', '#7f7f7f50', '#7f7f7f90'],
        autopct = lambda x: '{:.0f}'.format(x*counts.sum()/100),
    )

    ax[0,0].pie(
        counts.sum(axis=1), radius=1, wedgeprops=dict(width=0.3, edgecolor='w'),
        labels=['MSL1 peaks', 'no MSL1 peaks'],
        colors=['#1f77b4', '#7f7f7f']
        )

    ax[0,0].set_title("Downregulated genes in MSL1-gamma269")

    cs = [Line2D([0], [0], color='#1f77b490', lw=4),
                    Line2D([0], [0], color='#1f77b450', lw=4),
                    Line2D([0], [0], color='#7f7f7f90', lw=4),
                    Line2D([0], [0], color='#7f7f7f50', lw=4)
                ]
    ax[0,0].legend(cs, ['X','non-X','X', 'non-X'], bbox_to_anchor=(0.9, 0.4))


    ax[1, 0].pie(
        counts_up.flatten(), radius=1-0.4,
        wedgeprops=dict(width=0.4, edgecolor='w'),
        colors=['#1f77b490', '#1f77b450', '#7f7f7f50', '#7f7f7f90'],
        autopct = lambda x: '{:.0f}'.format(x*counts_up.sum()/100),
    )

    ax[1,0].pie(
        counts_up.sum(axis=1), radius=1, wedgeprops=dict(width=0.3, edgecolor='w'),
        labels=['MSL1 peaks', 'no MSL1 peaks'],
        colors=['#1f77b4', '#7f7f7f']
        )

    ax[1,0].set_title("Upregulated genes in MSL1-gamma269")

    cs = [Line2D([0], [0], color='#1f77b490', lw=4),
                    Line2D([0], [0], color='#1f77b450', lw=4),
                    Line2D([0], [0], color='#7f7f7f90', lw=4),
                    Line2D([0], [0], color='#7f7f7f50', lw=4)
                ]
    ax[1,0].legend(cs, ['X','non-X','X', 'non-X'], bbox_to_anchor=(0.9, 0.4))

    ######################################## Bar plots ########################################
    ####################################### down genes ########################################

    totx = sum(
        dict(MSL1[MSL1['gene_id'].isin(downgenes_withpeaks_X)]['annFeature'].value_counts()).values()
    )
    totnonx = sum(
        dict(MSL1[MSL1['gene_id'].isin(downgenes_withpeaks_nonX)]['annFeature'].value_counts()).values()
    )

    _c = {}
    for key in dict(MSL1[MSL1['gene_id'].isin(downgenes_withpeaks_X)]['annFeature'].value_counts()):
        _c[key] = [
            dict(MSL1[MSL1['gene_id'].isin(downgenes_withpeaks_X)]['annFeature'].value_counts())[key]/totx,
            dict(MSL1[MSL1['gene_id'].isin(downgenes_withpeaks_nonX)]['annFeature'].value_counts())[key]/totnonx
        ]
    _c

    chrom = (
        "X",
        "non-X"
    )

    width = 0.4
    bottom = np.zeros(2)

    for peaktype, _w in _c.items():
        _col=cdic[peaktype]
        ax[0,1].bar(chrom, _w, width, label=peaktype, bottom=bottom, color=_col)
        bottom += _w

    cs = []
    for peaktype in _c:
        cs.append(
            Line2D([0], [0], color=cdic[peaktype], lw=4)
        )

    ax[0,1].legend(cs, _c.keys())
    ax[0,1].set_ylabel("Fraction")
    ax[0,1].set_title("MSL1 peaks associated to downregulated genes")

    ####################################### up genes ########################################

    totx_up = sum(
        dict(MSL1[MSL1['gene_id'].isin(upgenes_withpeaks_X)]['annFeature'].value_counts()).values()
    )
    totnonx_up = sum(
        dict(MSL1[MSL1['gene_id'].isin(upgenes_withpeaks_nonX)]['annFeature'].value_counts()).values()
    )

    _c_up = {}
    for key in dict(MSL1[MSL1['gene_id'].isin(upgenes_withpeaks_X)]['annFeature'].value_counts()):
        _c_up[key] = [
            dict(MSL1[MSL1['gene_id'].isin(upgenes_withpeaks_X)]['annFeature'].value_counts())[key]/totx_up,
            dict(MSL1[MSL1['gene_id'].isin(upgenes_withpeaks_nonX)]['annFeature'].value_counts())[key]/totnonx_up
        ]
    _c_up

    width = 0.4
    bottom = np.zeros(2)

    for peaktype, _w in _c_up.items():
        _col=cdic[peaktype]
        ax[1,1].bar(chrom, _w, width, label=peaktype, bottom=bottom, color=_col)
        bottom += _w

    cs = []
    for peaktype in _c_up:
        cs.append(
            Line2D([0], [0], color=cdic[peaktype], lw=4)
        )

    ax[1,1].legend(cs, _c_up.keys())
    ax[1,1].set_ylabel("Fraction")
    ax[1,1].set_title("MSL1 peaks associated to upregulated genes")
    fig.savefig(ofile, dpi=300)

plotter(
    snakemake.input.finhits,
    snakemake.params.DE_MSL1,
    snakemake.params.padj_cutoff,
    snakemake.params.gene,
    snakemake.output[0]
)