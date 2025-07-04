include: 'rules/rna.smk'

rule all:
  input:
    'RNA/RNAseq.html',
    'RNA/upreg_deln.tsv',