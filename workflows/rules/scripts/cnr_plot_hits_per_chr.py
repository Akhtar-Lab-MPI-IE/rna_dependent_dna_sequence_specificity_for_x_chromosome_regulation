import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import gimmemotifs
from gimmemotifs.motif import read_motifs
from gimmemotifs.scanner import Scanner


pfmfile = snakemake.input.mot
chromsizes = snakemake.params.chromsizes
padj_cutoff = snakemake.params.cutoff
fimo = snakemake.input.fimo
opdf = snakemake.output.opdf

MOTORDER  = ['Cluster_9', 'Cluster_2', 'Cluster_4', 'Cluster_1' ,'Cluster_6', 'Cluster_3', 'Cluster_8', '4', '64', '63']
CHROMS = ['2L', '2R', '3L', '3R', 'X', 'Y']

a = pd.read_table(fimo)
a = a[a['p-value'] < padj_cutoff]
a = a[a['sequence_name'].isin(CHROMS)]

chromdic = {}
with open(chromsizes) as f:
    for line in f:
        chromdic[line.strip().split()[0]] = line.strip().split()[1]

fig, axs = plt.subplots(figsize=(12,6), constrained_layout=True)
axs.set_axis_off()

pfms = dict([(m.id, m) for m in read_motifs(pfmfile)])
motifs = [m for m in pfms.keys()]


s = Scanner()
s.set_motifs(pfmfile)

gs = fig.add_gridspec(7, len(motifs)-5)

fig, axs = plt.subplots(figsize=(12,6), constrained_layout=True)
axs.set_axis_off()

pfms = dict([(m.id, m) for m in read_motifs(pfmfile)])
motifs = [m for m in pfms.keys()]


s = Scanner()
s.set_motifs(pfmfile)

gs = fig.add_gridspec(7, len(motifs)-5)

for i, mot in enumerate(MOTORDER):
    if i <= 4:
        barax = fig.add_subplot(gs[0:2, i])
        _x = fig.add_subplot(gs[2:3, i])
    else:
        barax = fig.add_subplot(gs[4:6, i-5])
        _x = fig.add_subplot(gs[6:7, i-5])
    
    # Plot Motif
    pfms[mot].plot_logo(ax=_x, title=False, ylabel=False)
    _x.set_axis_off()
    
    
    # Plot bars
    counts = a[a['motif_id'] == mot].groupby('sequence_name').size()
    counts = pd.DataFrame(counts)
    counts.columns = ['counts']
    counts['chromsize'] = [int(chromdic[i])/1000000 for i in counts.index]
    counts['motfrac'] = counts['counts'] / counts['chromsize']
    counts
    sns.barplot(
        data=counts,
        x=counts.index,
        y=counts.motfrac,
        ax=barax
    )
    barax.set_xlabel("")
    barax.set_ylabel("Motifs / Mb")
    
    


plt.savefig(opdf, dpi=300)