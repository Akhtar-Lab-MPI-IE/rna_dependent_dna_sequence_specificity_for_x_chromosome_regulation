from gimmemotifs.motif import read_motifs
from gimmemotifs.scanner import Scanner
from gimmemotifs.fasta import Fasta
from gimmemotifs.plot import diff_plot
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
from matplotlib.colors import Normalize, to_hex
import seaborn as sns
import pandas as pd
from matplotlib.collections import PatchCollection
from matplotlib.colors import TwoSlopeNorm
import matplotlib.colors as mcolors
from matplotlib.legend_handler import HandlerPatch
import matplotlib.patches as mstatannotpatches
from matplotlib.ticker import MaxNLocator
from scipy.stats import fisher_exact, poisson_means_test
from statannotations.Annotator import Annotator

infiles = sorted(snakemake.input.fnas)
bgfile = snakemake.input.fnashuf
pfmfile = snakemake.input.motifs
padj_cutoff = snakemake.params.cutoff
chromsizes = snakemake.params.chromsizes
fimo = snakemake.input.fimo

catf_x = snakemake.output.catf_x
catf_a = snakemake.output.catf_a
ofile = snakemake.output.ofile
otab = snakemake.output.otab



chromsizes = '/data/repository/organisms/dm6_ensembl/genome_fasta/genome.chrom.sizes'
fimo = 'CutnRun/motifs/clusteredMotifs_fimo/fimo.tsv'
CHROMS = ['2L', '2R', '3L', '3R', '4', 'X']

cutoff = 0.9
minenr = 0
minfreq = 0
SIZE_NORMFAC = 1000
MINSC = 0
MAXSC = 0.2

# Combine fasta files (temporary)
for i in [i for i in infiles if '_X_' in i]:
    with open(catf_x, 'w') as f:
        with open(i, 'r') as red:
            for line in red:
                f.write(line)

for i in [i for i in infiles if '_autosome_' in i]:
    with open(catf_a, 'w') as f:
        with open(i, 'r') as red:
            for line in red:
                f.write(line)

# Manually defined, might change if the input motifs for clustering changes as well.
rename_dic = {
    'Cluster_6': 'Unknown - CAGCAACAACAACAACA',
    'Cluster_8': 'GAGA',
    'Cluster_4': 'M1BP',
    'Cluster_3': 'Unknown - GTGTGTGTGTGTG',
    'Cluster_2': 'Unknown - CACACACACACACAC',
    'Cluster_1': 'KLF10',
    'Cluster_9': 'DREF',
    '4': 'Unknown - nnnnnnTGGCAgCnC',
    '63': 'NTC Su(H)',
    '64': 'NTC sna',
}

ORDER = [
    'Cluster_9', 'Cluster_2', 'Cluster_4', 'Cluster_1', 'Cluster_6', 'Cluster_3', 'Cluster_8', '4', '64', '63'
]

pfms = dict([(m.id, m) for m in read_motifs(pfmfile)])
motifs = [m for m in pfms.keys()]
names = [fname.split('/')[-2] for fname in [catf_x, catf_a]]

s = Scanner()
s.set_motifs(pfmfile)
s.set_threshold(threshold=cutoff)

# Get background frequencies
nbg = float(len(Fasta(bgfile).seqs))
bgcounts = s.total_count(bgfile, nreport=1)

nx = float( len(Fasta(catf_x).seqs) )
xcounts = s.total_count(catf_x, nreport=1)

na = float( len(Fasta(catf_a).seqs) )
acounts = s.total_count(catf_a, nreport=1)

rawdf = pd.DataFrame(
    [xcounts, acounts, bgcounts]
)
rawdf.columns = motifs
rawdf.index = ['x', 'autosomes', 'bg']
rawdf['n'] = [nx, na, nbg]
rawdf

# Collect frequencies - fisher exact statistic, p-val

statres = []
for mot in motifs:
    xf = rawdf[mot].loc['x']
    xn = rawdf['n'].loc['x']
    
    af = rawdf[mot].loc['autosomes']
    an = rawdf['n'].loc['autosomes']
    
    _tabs = np.array(
        [
            (xf, xn - xf),
            (af, an - af)
        ]
    )
    
    _t = fisher_exact(_tabs)
    
    statres.append(
        [mot, xf/xn, af/an, _t[0], _t[1], _t[1] * len(motifs)]
    )

statres = pd.DataFrame(statres)
statres.columns = ['mot', 'X', 'Autosomes', 'stat', 'p-value', 'padj']
statres['mot'] =  pd.Categorical(statres['mot'], categories=ORDER, ordered=True)
statres = statres.sort_values('mot')
statres['mot'] = statres['mot'].map(rename_dic)
pdf = statres.melt(id_vars=['stat', 'p-value', 'padj', 'mot'])
statres.to_csv(otab, sep='\t', index=False)
a = pd.read_table(fimo)
a = a[a['p-value'] < padj_cutoff]
a = a[a['sequence_name'].isin(CHROMS)]
a
chromdic = {}
with open(chromsizes) as f:
    for line in f:
        chromdic[line.strip().split()[0]] = line.strip().split()[1]
chromdic
lensdic = {}
lensdic['X'] = int(chromdic['X'])
lensdic['autosomes'] = int(chromdic['2L']) + int(chromdic['2R']) + int(chromdic['3L']) + int(chromdic['3R']) + int(chromdic['4'])

countsres = []

for i, mot in enumerate(ORDER):
    # X
    xcounts = len(a[(a['motif_id'] == mot) & (a['sequence_name'] == 'X')])
    # Autosomes
    acounts = len(a[(a['motif_id'] == mot) & (a['sequence_name'] != 'X')])
    stat, pval = poisson_means_test(
        xcounts,
        lensdic['X'],
        acounts,
        lensdic['autosomes']
    )
    countsres.append(
        [mot, xcounts, xcounts/(lensdic['X']/1000000), acounts, acounts/(lensdic['autosomes']/1000000), stat, pval, pval*len(ORDER)]
    )
countdf = pd.DataFrame(countsres)
countdf.columns = ['mot', 'xcounts', 'xrel', 'acounts', 'arel', 'stat', 'p-val', 'padj']
countdf['mot'] =  pd.Categorical(countdf['mot'], categories=ORDER, ordered=True)
countdf = countdf.sort_values('mot')
countdf['mot_renamed'] = countdf['mot'].map(rename_dic)
# Create figure
fig = plt.figure(figsize=(14, 8), constrained_layout=True)  # Adjust width and height as needed
gs = GridSpec(nrows=2, ncols=10, figure=fig, height_ratios=[1, 1])
freqax = fig.add_subplot(gs[0, :])

sns.barplot(
    data=pdf,
    x='mot',
    y='value',
    hue='variable',
    ax=freqax
)
_pvals = []
_pairs = []
for i,r in statres.iterrows():
    if r['padj'] < 0.05:
        _pvals.append(r['padj'])
        _pairs.append(
            (
                (r['mot'], 'X'),(r['mot'], 'Autosomes')
            )
        )

annot = Annotator(
    freqax,
    pairs=_pairs,
    data=pdf,
    x="mot",
    y="value",
    hue="variable",
    plot="barplot",
)
annot.new_plot(ax=freqax, pairs=_pairs,
               data=pdf, x='mot', y='value', hue='variable')

(annot
 .configure(test=None)
 .set_pvalues(pvalues=_pvals)
 .annotate())
freqax.set_xlabel('')
freqax.set_ylabel('motifs/peak')
freqax.set_xticklabels(freqax.get_xticklabels(), rotation=90, ha='right')
freqax.legend(title='')

bpaxes = []
for i, mot in enumerate(ORDER):
    print(mot)
    ax = fig.add_subplot(gs[1, i])
    
    tdf = countdf[countdf['mot'] == mot]

    padj = tdf['padj'].values[0]
    plotdf = tdf.melt(id_vars='mot', value_vars=['xrel', 'arel'])
    
    sns.barplot(
        data=plotdf,
        x='variable',
        y='value',
        ax=ax
    )
    if padj < 0.05:
        _pairs = [('xrel', 'arel')]
        _pvals = [padj]
        
        annot = Annotator(
            ax,
            pairs=_pairs,
            data=plotdf,
            x="variable",
            y="value",
            plot="barplot",
        )
        
        annot.new_plot(ax=ax, pairs=_pairs,
               data=plotdf, x='variable', y='value')
        (annot
             .configure(test=None)
             .set_pvalues(pvalues=_pvals)
             .annotate())
        
    ax.set_xticklabels(['X', 'auto'])
    ax.set_xlabel(f"{mot}")
    if i == 0:
        ax.set_ylabel("Motifs/Mbp")
    else:
        ax.set_ylabel("")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, prune=None, nbins='auto'))
    bpaxes.append(ax)

fig.savefig(ofile, dpi=300)