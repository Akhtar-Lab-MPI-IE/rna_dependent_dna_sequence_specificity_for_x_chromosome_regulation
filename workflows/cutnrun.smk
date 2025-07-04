SAMPLES, = glob_wildcards("fqs/CNR/{sample}_R1.fastq.gz")
print(f"Samples = {SAMPLES}")

IP = ['CNR_MSL1']
peakCategories = [
  'MSL1_X_intergenic',
  'MSL1_X_promoter',
  'MSL1_X_gene',
  'MSL1_autosome_intergenic',
  'MSL1_autosome_promoter',
  'MSL1_autosome_gene',
]

chipDict = {
    "CNR_MSL1_rep1" : "CNR_IgG_rep1",
    "CNR_MSL1_rep2": "CNR_IgG_rep2"
}

include: 'rules/cutnrun_preprocessing.smk'
include: 'rules/cutnrun_peaks.smk'
include: 'rules/cutnrun_heatmap.smk'
include: 'rules/cutnrun_motifs.smk'
wildcard_constraints:
  sample = "[^.]+"


rule all:
  input:
    # preprocessing
    expand(
      'CutnRun/CutnRun_dedup/filtered_bam/{sample}.filtered.bam',
      sample=SAMPLES
    ),
    expand(
      'CutnRun/bws/{sample}.bw', sample=SAMPLES
    ),
    expand(
      'CutnRun/bws/{chip}.log2.bw', chip=chipDict.keys()
    ),
    # peaks
    expand(
      'CutnRun/peaks/{chip}.macs2_peaks.xls', chip = chipDict.keys()
    ),
    'CutnRun/final_peaks/CNR_MSL1.bed',
    'CutnRun/peakannotation/CNR_MSL1_finalhits.txt',
    'figures/cutnrun_peakannotation.pdf',
    'figures/cutnrun_rna.pdf',
    'figures/cutnrun_heatmap_MSL1.pdf',
    # Motif enrichments'
    expand('CutnRun/motifs/{peakcat}/meme', peakcat=peakCategories),
    expand('CutnRun/motifs/{peakcat}/homer', peakcat=peakCategories),
    'figures/MSL1peaks_motifs_perchr.pdf',
    'figures/MSL1peaks_motif_enrichment.pdf',
    'figures/MSL1peaks_motif_stats.pdf'