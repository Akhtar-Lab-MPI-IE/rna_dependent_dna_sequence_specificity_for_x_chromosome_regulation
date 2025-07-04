rule CNR_bed2fna:
  input:
    bed = 'CutnRun/motifs/{peakcat}/peaks.bed',
  output:
    fna = 'CutnRun/motifs/{peakcat}/peaks.fna',
  params:
    gfile = config['fna']
  shell:'''
  bedtools getfasta -fi {params.gfile} -bed {input.bed} > {output.fna}
  '''

rule CNR_dust:
  input:
    fna = 'CutnRun/motifs/{peakcat}/peaks.fna',
  output:
    fna = 'CutnRun/motifs/{peakcat}/peaks_dusted.fna',
  threads: 2
  conda: 'envs/meme.yml'
  shell:'''
  dust {input.fna} > {output.fna}
  '''

rule CNR_meme:
  input:
    fna = 'CutnRun/motifs/{peakcat}/peaks_dusted.fna',
  output:
    directory('CutnRun/motifs/{peakcat}/meme'),
  threads: 2
  conda: 'envs/meme.yml'
  shell:'''
  meme-chip -minw 5 -meme-searchsize 0 -meme-nmotifs 25 -o {output} {input.fna}
  '''

rule CNR_homer:
  input:
    fnaf = 'CutnRun/motifs/{peakcat}/peaks.fna',
  output:
    directory('CutnRun/motifs/{peakcat}/homer')
  params:
    gfile = config['fna']
  threads: 3
  conda: 'envs/homer.yml'
  shell:'''
  findMotifs.pl {input.fnaf} fasta {output} -fasta {params.gfile} -mask
  '''

rule CNR_shuffle:
  input:
    fnas = expand('CutnRun/motifs/{peakcat}/peaks.fna', peakcat = peakCategories)
  output:
    fnaout = temp('CutnRun/motifs/MSL1peaks.fna'),
    fnashuf = 'CutnRun/motifs/MSL1peaks_shuffled.fna'
  threads: 2
  conda: 'envs/meme.yml'
  shell:'''
  cat {input.fnas} > {output.fnaout}
  fasta-shuffle-letters -seed 12345 {output.fnaout} > {output.fnashuf}
  '''

rule CNR_collate_motifs:
  input:
    expand('CutnRun/motifs/{peakcat}/meme', peakcat=peakCategories),
    expand('CutnRun/motifs/{peakcat}/homer', peakcat=peakCategories)
  output:
    mot = 'CutnRun/motifs/enrichedmotifs.meme'
  params:
    hdir = 'CutnRun/motifs/',
    extra_mot = config['motifs'],
    cutoff = 1e-15
  threads: 2
  script: 'scripts/cnr_collate_motifs.py'

rule CNR_cluster_motifs:
  input:
    mot = 'CutnRun/motifs/enrichedmotifs.meme'
  output:
    mot = 'CutnRun/motifs/clusteredMotifs/enriched_motifs_consensus_motifs.meme'
  conda: 'envs/gimmemotifs.yml'
  threads: 10
  shell:'''
  TOBIAS ClusterMotifs -t 0.4 -a meme -o CutnRun/motifs/clusteredMotifs \
    -p enriched_motifs -m {input.mot}
  '''

rule CNR_fimo_clustermotifs:
  input:
    mot = 'CutnRun/motifs/clusteredMotifs/enriched_motifs_consensus_motifs.meme'
  output:
    fimof = 'CutnRun/motifs/clusteredMotifs_fimo/fimo.tsv'
  params:
    fna = config['fna'],
    oc = 'CutnRun/motifs/clusteredMotifs_fimo'
  conda: 'envs/meme.yml'
  threads: 2
  shell:'''
  fimo --oc {params.oc} {input.mot} {params.fna}
  '''

rule CNR_clusmot_to_mot:
  input:
    mot = 'CutnRun/motifs/clusteredMotifs/enriched_motifs_consensus_motifs.meme'
  output:
    mot = 'CutnRun/motifs/clusteredMotifs/enriched_motifs_consensus_motifs.mot'
  threads: 2
  script: 'scripts/cnr_clusmot_to_gimmemotifs.py'

rule CNR_fimo_plot:
  input:
    mot = 'CutnRun/motifs/clusteredMotifs/enriched_motifs_consensus_motifs.mot',
    fimo = 'CutnRun/motifs/clusteredMotifs_fimo/fimo.tsv'
  output:
    opdf = 'figures/MSL1peaks_motifs_perchr.pdf'
  params:
    chromsizes = config['chromsizes'],
    cutoff = 1e-5
  conda: 'envs/gimmemotifs.yml'
  script: 'scripts/cnr_plot_hits_per_chr.py'

rule CNR_motenrich:
  input:
    fnas = expand('CutnRun/motifs/{peakcat}/peaks_dusted.fna', peakcat = peakCategories),
    fnashuf = 'CutnRun/motifs/MSL1peaks_shuffled.fna',
    motifs = 'CutnRun/motifs/clusteredMotifs/enriched_motifs_consensus_motifs.mot'
  output:
    ofile = 'figures/MSL1peaks_motif_enrichment.pdf'
  conda: 'envs/gimmemotifs.yml'
  script: 'scripts/cnr_motifenrichment.py'

rule CNR_mots_xvsauto:
  input:
    fnas = expand('CutnRun/motifs/{peakcat}/peaks_dusted.fna', peakcat = peakCategories),
    fnashuf = 'CutnRun/motifs/MSL1peaks_shuffled.fna',
    motifs = 'CutnRun/motifs/clusteredMotifs/enriched_motifs_consensus_motifs.mot',
    fimo = 'CutnRun/motifs/clusteredMotifs_fimo/fimo.tsv'
  output:
    catf_x = temp('CutnRun/motifs/MSL1peaks_X.fna'),
    catf_a = temp('CutnRun/motifs/MSL1peaks_auto.fna'),
    ofile = 'figures/MSL1peaks_motif_stats.pdf',
    otab = 'CutnRun/motifs/enrichment_stats.tsv'
  params:
    chromsizes = config['chromsizes'],
    cutoff = 1e-5
  conda: 'envs/gimmemotifs.yml'
  script: 'scripts/cnr_mots_xvsauto.py'
