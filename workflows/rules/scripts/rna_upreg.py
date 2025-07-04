import pandas as pd

delc = snakemake.params.delc
dell1 = snakemake.params.dell1
dell2 = snakemake.params.dell2
deln = snakemake.params.deln
msl1null = snakemake.params.msl1null
padj_cutoff = snakemake.params.padj_cutoff

delc_out = snakemake.output.delc
dell1_out = snakemake.output.dell1
dell2_out = snakemake.output.dell2
deln_out = snakemake.output.deln

# Get upregulated genes 
msl1null = pd.read_table(msl1null, sep='\t')
msl1null = msl1null[(msl1null['padj'] < padj_cutoff) & (msl1null['log2FoldChange'] > 0)]
null_up_genes = set(msl1null.index)
null_up_genes

# delc
delc = pd.read_table(delc, sep='\t')
delc = delc[(delc['padj'] < padj_cutoff) & (delc['log2FoldChange'] > 0)]
delc = delc[~delc.index.isin(null_up_genes)]
delc.to_csv(delc_out, sep='\t')
# dell1
dell1 = pd.read_table(dell1, sep='\t')
dell1 = dell1[(dell1['padj'] < padj_cutoff) & (dell1['log2FoldChange'] > 0)]
dell1 = dell1[~dell1.index.isin(null_up_genes)]
dell1.to_csv(dell1_out, sep='\t')
# dell2
dell2 = pd.read_table(dell2, sep='\t')
dell2 = dell2[(dell2['padj'] < padj_cutoff) & (dell2['log2FoldChange'] > 0)]
dell2 = dell2[~dell2.index.isin(null_up_genes)]
dell2.to_csv(dell2_out, sep='\t')
# deln
deln = pd.read_table(deln, sep='\t')
deln = deln[(deln['padj'] < padj_cutoff) & (deln['log2FoldChange'] > 0)]
deln = deln[~deln.index.isin(null_up_genes)]
deln.to_csv(deln_out, sep='\t')