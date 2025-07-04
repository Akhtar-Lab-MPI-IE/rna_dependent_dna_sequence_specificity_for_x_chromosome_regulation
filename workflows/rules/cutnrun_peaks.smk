rule CNR_macs2_igg:
  input:
    bams = expand('CutnRun/sieve/{sample}.bam', sample=SAMPLES),
  output:
    'CutnRun/peaks/{chip}.macs2_peaks.xls'
  params:
    pref = lambda wildcards: 'CutnRun/peaks/' + wildcards.chip + '.macs2',
    chip = lambda wildcards: 'CutnRun/sieve/' + wildcards.chip + ".bam",
    ctrl = lambda wildcards: 'CutnRun/sieve/' + chipDict[wildcards.chip] + ".bam"
  conda:
    'envs/macs.yml'
  shell:'''
  macs3 callpeak -t {params.chip} \
    -p 1e-3 -f BAM --keep-dup all -n {params.pref} \
    --cutoff-analysis
  '''

rule CNR_filter_peaks:
  input:
    'CutnRun/peaks/{chip}.macs2_peaks.xls'
  output:
    canchr = temp('CutnRun/final_peaks/{chip}.can.bed'),
    grepf = temp('CutnRun/final_peaks/{chip}_can_chrs.txt')
  params:
    peak = lambda wildcards: 'CutnRun/peaks/' + wildcards.chip + '.macs2_peaks.narrowPeak',
    can_chr = config['canchr']
  shell:'''
  awk '{{print "^" $0}}' {params.can_chr} > {output.grepf}
  grep -f {output.grepf} {params.peak} > {output.canchr}
  '''

rule CNR_intersect_peaks:
  input:
    MSL1_r1 = temp("CutnRun/final_peaks/CNR_MSL1_rep1.can.bed"),
    MSL1_r2 = temp("CutnRun/final_peaks/CNR_MSL1_rep2.can.bed")
  output:
    MSL1 = 'CutnRun/final_peaks/CNR_MSL1.bed'
  shell:'''
  bedtools intersect -a {input.MSL1_r1} -b {input.MSL1_r2} > {output.MSL1}
  '''

rule CNR_annotate_peaks:
  input:
    peakset = 'CutnRun/final_peaks/CNR_MSL1.bed',
  output:
    finhits = 'CutnRun/peakannotation/CNR_MSL1_finalhits.txt',
    finhitsbed = temp('CutnRun/peakannotation/CNR_MSL1_finalhits.bed'),
    allhitsbed = temp('CutnRun/peakannotation/CNR_MSL1_allhits.bed'),
    allhitstxt = temp('CutnRun/peakannotation/CNR_MSL1_allhits.txt'),
    json = temp('CutnRun/peakannotation/CNR_MSL1.json'),
    sumr = temp('CutnRun/peakannotation/CNR_MSL1_summary.pdf')
  conda: 'envs/uropa.yml'
  params:
    gtf = 'RNA/RNAseq/Annotation/genes.filtered.gtf',
    odir = 'CutnRun/peakannotation/',
  threads: 10
  shell:'''
  uropa -b {input.peakset} \
    -g {params.gtf} \
    --summary --feature transcript,three_prime_utr,five_prime_utr,exon \
    -p CNR_MSL1 \
    --distance 10000 5000 \
    --internals 1 -o {params.odir} \
    -t {threads} --show-attributes gene_id transcript_id gene_name transcript_biotype
  '''

rule CNR_plot_peakanno:
  input:
    MSL1 = 'CutnRun/peakannotation/CNR_MSL1_finalhits.txt'
  output:
    'figures/cutnrun_peakannotation.pdf'
  params:
    gene = config['gtf']
  threads: 2
  script:
    'scripts/cnr_plot_peakanno.py'

rule CNR_plot_peakgeneconnection:
  input:
    finhits = 'CutnRun/peakannotation/CNR_MSL1_finalhits.txt'
  output:
    'figures/cutnrun_rna.pdf'
  params:
    DE_MSL1 = 'RNA/RNA_WT_vs_msl1null_shrunk.tsv',
    padj_cutoff = config['padj_cutoff'],
    gene = config['gtf']
  threads: 2
  script:
    'scripts/cnr_plot_peakannogene.py'

rule CNR_splitpeaks:
  input:
    peakset = 'CutnRun/peakannotation/CNR_MSL1_finalhits.txt'
  output:
    ofx1 = 'CutnRun/motifs/MSL1_X_intergenic/peaks.bed',
    ofx2 = 'CutnRun/motifs/MSL1_X_promoter/peaks.bed',
    ofx3 = 'CutnRun/motifs/MSL1_X_gene/peaks.bed',
    ofnonx1 = 'CutnRun/motifs/MSL1_autosome_intergenic/peaks.bed',
    ofnonx2 = 'CutnRun/motifs/MSL1_autosome_promoter/peaks.bed',
    ofnonx3 = 'CutnRun/motifs/MSL1_autosome_gene/peaks.bed',
  threads: 2
  script:
    'scripts/cnr_splitpeaks.py'

