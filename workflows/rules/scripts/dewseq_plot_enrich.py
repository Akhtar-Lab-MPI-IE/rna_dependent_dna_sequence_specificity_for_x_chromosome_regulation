import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

cmap = {
    'roX1': "#ff7f0e",
    'roX2': "#1f77b4",
    'snoRNA:185': "#8c564b",
    'snoRNA:U29:54Ed': "#9467bd",
    'Uhg1': "#d62728",
    'non-significant': '#7f7f7f'
}

def parse_dewseqtable(df, FDRcut):
    _t = []
    for i,r in df.iterrows():
        if r['pSlidingWindows.adj'] < FDRcut:
            _t.append(r['gene_name'])
        else:
            _t.append('non-significant')
    df['sigs'] = _t
    return(df)

FDRcut = snakemake.params.FDRcut
C8 = pd.read_table(snakemake.params.C8)
C8 = parse_dewseqtable(C8, FDRcut)
S2 = pd.read_table(snakemake.params.S2)
S2 = parse_dewseqtable(S2, FDRcut)
MLE = pd.read_table(snakemake.params.MLE)
MLE = parse_dewseqtable(MLE, FDRcut)
MSL2 = pd.read_table(snakemake.params.MSL2)
MSL2 = parse_dewseqtable(MSL2, FDRcut)
MLEf = pd.read_table(snakemake.params.MLEf)
MLEf = parse_dewseqtable(MLEf, FDRcut)
MSL2f = pd.read_table(snakemake.params.MSL2f)
MSL2f = parse_dewseqtable(MSL2f, FDRcut)

fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(15,18), sharey=True)

ax=axs[0,0]
sns.scatterplot(
    data=C8,
    x='log2FoldChange',
    y=-np.log10(C8['pSlidingWindows.adj']),
    hue=C8['sigs'],
    edgecolor='none',
    palette=cmap,
    hue_order=cmap.keys(),
    ax=ax
)
ax.set_title("C8")
ax.set_ylabel("-log10(Padj)")
ax.legend(title="")

ax=axs[0,1]
sns.scatterplot(
    data=S2,
    x='log2FoldChange',
    y=-np.log10(S2['pSlidingWindows.adj']),
    hue=S2['sigs'],
    edgecolor='none',
    palette=cmap,
    hue_order=cmap.keys(),
    ax=ax
)
ax.set_title("S2")
ax.get_legend().set_visible(False)

ax=axs[1,0]
sns.scatterplot(
    data=MLE,
    x='log2FoldChange',
    y=-np.log10(MLE['pSlidingWindows.adj']),
    hue=MLE['sigs'],
    edgecolor='none',
    palette=cmap,
    hue_order=cmap.keys(),
    ax=ax
)
ax.set_title("MLE")
ax.get_legend().set_visible(False)
ax.set_ylabel("-log10(Padj)")

ax=axs[1,1]
sns.scatterplot(
    data=MSL2,
    x='log2FoldChange',
    y=-np.log10(MSL2['pSlidingWindows.adj']),
    hue=MSL2['sigs'],
    edgecolor='none',
    palette=cmap,
    hue_order=cmap.keys(),
    ax=ax
)
ax.set_title("MSL2")
ax.get_legend().set_visible(False)

ax=axs[2,0]
sns.scatterplot(
    data=MLEf,
    x='log2FoldChange',
    y=-np.log10(MLEf['pSlidingWindows.adj']),
    hue=MLEf['sigs'],
    edgecolor='none',
    palette=cmap,
    hue_order=cmap.keys(),
    ax=ax
)
ax.set_title("MLE - male vs female")
ax.get_legend().set_visible(False)
ax.set_ylabel("-log10(Padj)")

ax=axs[2,1]
sns.scatterplot(
    data=MSL2f,
    x='log2FoldChange',
    y=-np.log10(MSL2f['pSlidingWindows.adj']),
    hue=MSL2f['sigs'],
    edgecolor='none',
    palette=cmap,
    hue_order=cmap.keys(),
    ax=ax
)
ax.set_title("MSL2 - male vs female")
ax.get_legend().set_visible(False)

fig.savefig(snakemake.output.opng, dpi=300)