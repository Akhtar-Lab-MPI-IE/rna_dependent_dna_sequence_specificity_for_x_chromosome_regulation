import glob

rule RNA_snakepipes:
  localrule: True
  output:
    'RNA/RNAseq/mRNAseq_snakePipes.done'
  threads: 1
  conda:
    config['snakePipesEnv']
  shell:'''
  mRNAseq -i fqs/RNA \
  --fastqc --trim -o RNA/RNAseq \
  dm6
  '''

rule RNA_DE:
  input:
    'RNA/RNAseq/mRNAseq_snakePipes.done'
  output:
    'RNA/RNAseq.html'
  threads: 2
  conda: 'envs/deseq.yml'
  params:
    counts = 'RNA/RNAseq/featureCounts/counts.tsv',
    t2g = 'RNA/RNAseq/Annotation/genes.filtered.t2g',
    odir = 'RNA/',
    geneinfo = config['geneinfo'],
    padj_cutoff = config['padj_cutoff'],
    pcaout = 'RNA/RNAseq_PCA.pdf',
    upsetout = 'RNA/RNAseq_upset.pdf'
  script:
    'scripts/RNA_DE.Rmd'

rule RNA_upreg:
  input:
    'RNA/RNAseq.html'
  output:
    delc = 'RNA/upreg_delc.tsv',
    dell1 = 'RNA/upreg_dell1.tsv',
    dell2 = 'RNA/upreg_dell2.tsv',
    deln = 'RNA/upreg_deln.tsv',
  threads: 2
  params:
    delc = 'RNA/RNA_WT_vs_delc_shrunk.tsv',
    dell1 = 'RNA/RNA_WT_vs_dell1_shrunk.tsv',
    dell2 = 'RNA/RNA_WT_vs_dell2_shrunk.tsv',
    deln = 'RNA/RNA_WT_vs_deln_shrunk.tsv',
    msl1null = 'RNA/RNA_WT_vs_delc_shrunk.tsv',
    padj_cutoff = config['padj_cutoff']
  script:
    'scripts/rna_upreg.py'