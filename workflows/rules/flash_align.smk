rule FLASH_fqc:
  input:
    R1 = 'fqs/FLASH/{sample}_R1.fastq.gz',
    R2 = 'fqs/FLASH/{sample}_R2.fastq.gz'
  output:
    'FLASH/fqc/{sample}_R1_fastqc.zip'
  params:
    odir = 'FLASH/fqc'
  conda: 'envs/flash.yml'
  threads: 5
  shell:'''
  fastqc -t {threads} -o {params.odir} {input.R1} {input.R2}
  '''


rule FLASH_star:
  input:
    R1 = 'fqs/FLASH/{sample}_R1.fastq.gz',
    R2 = 'fqs/FLASH/{sample}_R2.fastq.gz',
    fqc = 'FLASH/fqc/{sample}_R1_fastqc.zip'
  output:
    bam = 'FLASH/STAR/{sample}.Aligned.sortedByCoord.out.bam'
  params:
    gtf = config['gtf'],
    gdir = config['starindex'],
    sample = lambda wildcards: wildcards.sample
  conda: 'envs/flash.yml'
  threads: 20
  shell:'''
  STAR --runThreadN {threads} \
    --limitBAMsortRAM 128000000000 \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --sjdbGTFfile {params.gtf} \
    --genomeDir {params.gdir} \
    --outFileNamePrefix FLASH/STAR/{params.sample}. \
    --outFilterScoreMinOverLread 0.3 \
    --outFilterMatchNminOverLread 0.3 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --outFilterMismatchNmax 999 \
    --scoreDelOpen -1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --readFilesIn {input.R1} {input.R2} \
    --outFilterMultimapNmax 1 \
    --alignEndsType Extend5pOfRead1
  '''

rule FLASH_starix:
  input:
    bam = 'FLASH/STAR/{sample}.Aligned.sortedByCoord.out.bam',
  output:
    bai = 'FLASH/STAR/{sample}.Aligned.sortedByCoord.out.bam.bai'
  threads: 5
  conda: 'envs/flash.yml'
  shell:'''
  samtools index -@ {threads} {input.bam}
  '''

rule FLASH_dedup:
  input:
    bam = 'FLASH/STAR/{sample}.Aligned.sortedByCoord.out.bam',
    bai = 'FLASH/STAR/{sample}.Aligned.sortedByCoord.out.bam.bai'
  output:
    fbam = 'FLASH/dedup/{sample}.full.bam',
    bam = 'FLASH/dedup/{sample}.bam',
    bai = 'FLASH/dedup/{sample}.bam.bai',
    log = 'FLASH/dedup/{sample}.dedup.log'
  threads: 5
  conda: 'envs/flash.yml'
  shell:'''
  umi_tools dedup --stdin={input.bam} --log={output.log} > {output.fbam}
  samtools index -@ {threads} {output.fbam}
  samtools view -b {output.fbam} 2L 2R 3L 3R 4 X Y > {output.bam}
  samtools index -@ {threads} {output.bam}
  '''

rule FLASH_peaks:
  input:
    bams = expand('FLASH/dedup/{sample}.bam', sample=SAMPLES)
  output:
    cl = 'FLASH/peaks/{ipsample}.CL.bed',
    peaks = 'FLASH/peaks/{ipsample}.bed',
  params:
    fna = config['fna'],
    sbam = lambda wildcards: 'FLASH/dedup/' + wildcards.ipsample + '.bam',
    cbam = lambda wildcards: 'FLASH/dedup/' + CTRLDIC[wildcards.ipsample] + '.bam'
  conda: 'envs/flash.yml'
  threads: 15
  shell:'''
  pureclip -i {params.sbam} -bai {params.sbam}.bai \
    -ibam {params.cbam} -ibai {params.cbam}.bai \
    -g {params.fna} -nt {threads} -o {output.cl} -or {output.peaks}
  '''

rule FLASH_multiqc:
  input:
    bams = expand('FLASH/dedup/{sample}.full.bam', sample=SAMPLES)
  output:
    'FLASH/multiqc/multiqc_report.html'
  params:
    odir = 'FLASH/multiqc'
  conda: 'envs/flash.yml'
  shell:'''
  multiqc -o {params.odir} FLASH
  '''

rule FLASH_annotate_peaks:
  input:
    peakset = 'FLASH/peaks/{ipsample}.CL.bed'
  output:
    finhits = 'FLASH/peakannotation/{ipsample}_finalhits.txt',
    finhitsbed = temp('FLASH/peakannotation/{ipsample}_finalhits.bed'),
    allhitsbed = temp('FLASH/peakannotation/{ipsample}_allhits.bed'),
    allhitstxt = temp('FLASH/peakannotation/{ipsample}_allhits.txt'),
    json = temp('FLASH/peakannotation/{ipsample}.json'),
    sumr = temp('FLASH/peakannotation/{ipsample}_summary.pdf')
  conda: 'envs/uropa.yml'
  params:
    gtf = 'RNA/RNAseq/Annotation/genes.filtered.gtf',
    odir = 'FLASH/peakannotation/',
    prefix = lambda wildcards: wildcards.ipsample
  threads: 10
  shell:'''
  uropa -b {input.peakset} \
    -g {params.gtf} \
    --summary --feature transcript,three_prime_utr,five_prime_utr,exon \
    -p {params.prefix} \
    --distance 10000 5000 \
    --internals 1 -o {params.odir} \
    -t {threads} --show-attributes gene_id transcript_id gene_name transcript_biotype
  '''

rule FLASH_plot_clanno:
  input:
    expand('FLASH/peakannotation/{ipsample}_finalhits.txt', ipsample = CTRLDIC.keys())
  output:
    opng = 'figures/FLASH_annotation.pdf'
  params:
    annodir = 'FLASH/peakannotation'
  script: 'scripts/flash_plot_clanno.py'

rule FLASH_bamCompare:
  input:
    b = expand('FLASH/peakannotation/{ipsample}_finalhits.txt', ipsample = CTRLDIC.keys()),
  output:
    bw = 'FLASH/bw/{ipsample}.log2.bw'
  params:
    sample = lambda wildcards: 'FLASH/dedup/' + wildcards.ipsample + '.bam',
    igg = lambda wildcards: 'FLASH/dedup/' + CTRLDIC[wildcards.ipsample] + '.bam'
  threads: 5
  shell:'''
  bamCompare -p {threads} -b1 {params.sample} -b2 {params.igg} -bs 1 -o {output.bw}
  '''

rule FLASH_bamCompare_mvf:
  input:
    b = expand('FLASH/peakannotation/{ipsample}_finalhits.txt', ipsample = CTRLDIC.keys()),
  output:
    bw = 'FLASH/bw_mf/{malef}.mvsf.log2.bw'
  params:
    sample = lambda wildcards: 'FLASH/dedup/' + wildcards.malef + '.bam',
    igg = lambda wildcards: 'FLASH/dedup/' + CTRLDIC[wildcards.malef] + '.bam'
  threads: 5
  shell:'''
  bamCompare -p {threads} -b1 {params.sample} -b2 {params.igg} -bs 1 -o {output.bw}
  '''

rule FLASH_bamCoverage:
  input:
    'FLASH/dedup/{sample}.bam'
  output:
    'FLASH/bw/{sample}.RPKM.bw'
  threads: 10
  resources:
    mem_mb = 1000
  shell:'''
  bamCoverage -b {input} -p {threads} --normalizeUsing RPKM -o {output}
  '''

rule FLASH_mbwsumRPKM:
  input:
    expand('FLASH/bw/{sample}.RPKM.bw', sample=SAMPLES)
  output:
    mat = 'FLASH/plotcorr/{CORR}_matrixRPKM.npz'
  params:
    bws = lambda wildcards: expand('FLASH/bw/{CORRSAMPLE}.RPKM.bw', CORRSAMPLE=CORDIC[wildcards.CORR]),
    labels = lambda wildcards: CORDIC[wildcards.CORR]
  threads: 10
  resources:
    mem_mb = 1000
  shell:'''
  multiBigwigSummary bins -p {threads} --bw {params.bws} -o {output.mat} --labels {params.labels}
  '''

rule FLASH_plotcorrRPKM:
  input:
    'FLASH/plotcorr/{CORR}_matrixRPKM.npz'
  output:
    'FLASH/plotcorr/{CORR}_plotcorrRPKM.pdf'
  threads: 2
  shell:'''
  plotCorrelation -in {input} -o {output} -c spearman -p heatmap --removeOutliers --zMin 0.5 --zMax 1
  '''
