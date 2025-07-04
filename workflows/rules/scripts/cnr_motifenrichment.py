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
import matplotlib.patches as mpatches

infiles = snakemake.input.fnas
bgfile = snakemake.input.fnashuf
ofile = snakemake.output.ofile
pfmfile = snakemake.input.motifs
cutoff = 0.9
minenr = 0
minfreq = 0
SIZE_NORMFAC = 1000
MINSC = 0 # fdf['freq'].min()
MAXSC = 0.2 # fdf['freq'].max()

# Manually defined, might change if the input motifs for clustering changes as well.
rename_dic = {
    'Cluster_6': 'Unknown',
    'Cluster_8': 'GAGA',
    'Cluster_4': 'M1BP',
    'Cluster_3': 'Unknown',
    'Cluster_2': 'Unknown',
    'Cluster_1': 'KLF10',
    'Cluster_9': 'DREF',
    '4': 'Unknown',
    '63': 'NTC Su(H)',
    '64': 'NTC sna',
}

pfms = dict([(m.id, m) for m in read_motifs(pfmfile)])
motifs = [m for m in pfms.keys()]
names = [fname.split('/')[-2] for fname in infiles]

s = Scanner()
s.set_motifs(pfmfile)
s.set_threshold(threshold=cutoff)

# Get background frequencies
nbg = float(len(Fasta(bgfile).seqs))

bgcounts = s.total_count(bgfile, nreport=1)
bgfreq = [(c + 0.01) / nbg for c in bgcounts]

# Get frequences in input files
freq = {}
counts = {}
for fname in infiles:
    mcounts = s.total_count(fname, nreport=1)
    n = float(len(Fasta(fname).seqs))
    counts[fname] = mcounts
    freq[fname] = [(c + 0.01) / n for c in mcounts]


freqdf = pd.DataFrame(freq)
freqdf['background'] = bgfreq
freqdf.index = motifs
# Order df according to frequencies, separate true hits and NTC's though.
NTC = []
TMS = []
for i in freqdf.index:
    if 'NTC' in rename_dic[i]:
        NTC.append(i)
    else:
        TMS.append(i)
ordered_tms = list(freqdf.loc[TMS].max(axis=1).sort_values(ascending=True).index)
ixOrder = list(freqdf.loc[NTC].max(axis=1).sort_values(ascending=True).index) + ordered_tms
freqdf = freqdf.loc[ixOrder]
freqdf.columns = [i.split('/')[-2].replace("MSL1_", "").replace("_", " ") if '/' in i else 'background' for i in list(freqdf.columns)]
freqdf = freqdf[['X promoter', 'X gene', 'X intergenic', 'autosome promoter', 'autosome gene', 'autosome intergenic', 'background']]
enrdf = np.log2(freqdf.div(freqdf['background'], axis=0))
freqdf['mot'] = freqdf.index
fdf = freqdf.melt(id_vars='mot')
fdf.columns = ['mot', 'peaks', 'freq']
enrdf['mot'] = enrdf.index
rdf = enrdf.melt(id_vars='mot')
rdf.columns = ['mot', 'peaks', 'enrichment']
fdf = fdf.merge(rdf, left_on=['mot', 'peaks'], right_on=['mot', 'peaks'])
fdf['freqnorm'] = (fdf['freq'] - MINSC)/(MAXSC - MINSC)*SIZE_NORMFAC


fig, axs = plt.subplots(figsize=(16,12), constrained_layout=True)
axs.set_axis_off()
gs = fig.add_gridspec(len(motifs),10)
ax_hm = fig.add_subplot(gs[0:, 2:10], frameon=False)

norm = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)
cmap = mcolors.LinearSegmentedColormap.from_list(
    "", ["#1f77b4", "gray", "#d62728"]
)
im = ax_hm.scatter(
    fdf['peaks'], fdf['mot'], s=fdf['freqnorm'], cmap=cmap, c=fdf['enrichment'], norm=norm
)
ax_hm.set_xlim(-0.51, 6.51)
ax_hm.set_ylim(-0.51,len(motifs) - 0.49)
ax_hm.set_xticks(np.arange(len(freqdf.columns)) - 0.5, minor=True)
ax_hm.set_yticks(np.arange(len(freqdf.index)+1)-0.5, minor=True)
ax_hm.set_yticklabels([rename_dic[i] for i in freqdf.index])
ax_hm.set_xticklabels(ax_hm.get_xticklabels(), rotation=90)
ax_hm.grid(which='minor')
cb = plt.colorbar(im, ax=ax_hm)
cb.set_label("Log2 enrichment over background")

p1v = (0.01-MINSC)/(MAXSC - MINSC)*SIZE_NORMFAC
p20v = (0.2-MINSC)/(MAXSC - MINSC)*SIZE_NORMFAC
p40v = (0.4-MINSC)/(MAXSC - MINSC)*SIZE_NORMFAC


p1 = plt.scatter([],[], s=p1v, marker='o', color='#555555')
p20 = plt.scatter([],[], s=p20v, marker='o', color='#555555')
p40 = plt.scatter([],[], s=p40v, marker='o', color='#555555')

ax_hm.legend((p1, p20, p40),
    ('1%', '20%', '40%'),
    scatterpoints=1,
    bbox_to_anchor=[1.25, 0.6],
    ncol=1,
    labelspacing = 3,
    title='% of peaks with motif',
    handletextpad=1.5,
    borderpad=2,
    frameon=False
)

# Plot motifs
for i,mot in enumerate(list(freqdf.index)[::-1]):
    _x = fig.add_subplot(gs[i, 0:2])
    pfms[mot].plot_logo(ax=_x, title=False, ylabel=False)
    _x.set_axis_off()

plt.savefig(ofile, dpi=300)



fig, axs = plt.subplots(figsize=(16,12), constrained_layout=True)
axs.set_axis_off()
gs = fig.add_gridspec(len(motifs),10)
ax_hm = fig.add_subplot(gs[0:, 2:10], frameon=False)

norm = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)
cmap = mcolors.LinearSegmentedColormap.from_list(
    "", ["#d62728", "gray", "#1f77b4"]
)
im = ax_hm.scatter(
    fdf['peaks'], fdf['mot'], s=fdf['freqnorm'], cmap=cmap, c=fdf['enrichment'], norm=norm
)
ax_hm.set_xlim(-0.51, 6.51)
ax_hm.set_ylim(-0.51,len(motifs) - 0.49)
ax_hm.set_xticks(np.arange(len(freqdf.columns)) - 0.5, minor=True)
ax_hm.set_yticks(np.arange(len(freqdf.index)+1)-0.5, minor=True)
ax_hm.set_yticklabels([rename_dic[i] for i in freqdf.index])
ax_hm.set_xticklabels(ax_hm.get_xticklabels(), rotation=90)
ax_hm.grid(which='minor')
cb = plt.colorbar(im, ax=ax_hm)
cb.set_label("Log2 enrichment over background")

p1v = (0.01-MINSC)/(MAXSC - MINSC)*SIZE_NORMFAC
p20v = (0.2-MINSC)/(MAXSC - MINSC)*SIZE_NORMFAC
p40v = (0.4-MINSC)/(MAXSC - MINSC)*SIZE_NORMFAC


p1 = plt.scatter([],[], s=p1v, marker='o', color='#555555')
p20 = plt.scatter([],[], s=p20v, marker='o', color='#555555')
p40 = plt.scatter([],[], s=p40v, marker='o', color='#555555')

ax_hm.legend((p1, p20, p40),
    ('1%', '20%', '40%'),
    scatterpoints=1,
    bbox_to_anchor=[1.25, 0.6],
    ncol=1,
    labelspacing = 3,
    title='% of peaks with motif',
    handletextpad=1.5,
    borderpad=2,
    frameon=False
)

# Plot motifs
for i,mot in enumerate(list(freqdf.index)[::-1]):
    _x = fig.add_subplot(gs[i, 0:2])
    pfms[mot].plot_logo(ax=_x, title=False, ylabel=False)
    _x.set_axis_off()

plt.savefig(ofile.replace('.pdf', '_inverse.pdf'), dpi=300)