rule CNR_snakepipes:
  localrule: True
  output:
    snakerun = 'CutnRun/CutnRun_dedup/DNAmapping_snakePipes.done',
    bams = expand('CutnRun/CutnRun_dedup/filtered_bam/{sample}.filtered.bam', sample=SAMPLES)
  threads: 1
  conda:
    config['snakePipesEnv']
  shell:'''
  DNAmapping -i fqs/CNR \
    -o CutnRun/CutnRun_dedup \
    --cutntag --dedup --fastqc --trim dm6
  '''

rule CNR_sieve:
  input:
    snakerun = 'CutnRun/CutnRun_dedup/DNAmapping_snakePipes.done'
  output:
    bam = 'CutnRun/sieve/{sample}.bam',
    bai = 'CutnRun/sieve/{sample}.bam.bai',
    qc = 'CutnRun/sieve/{sample}.txt'
  params:
    rar = config['rar'],
    bam = lambda wildcards: f"CutnRun/CutnRun_dedup/filtered_bam/{wildcards.sample}.filtered.bam"
  threads: 10
  shell:'''
  alignmentSieve --bam {params.bam} --outFile {output.bam} \
    --filterMetrics {output.qc} -p {threads} --blackListFileName {params.rar} \
    --maxFragmentLength 1000 --minFragmentLength 0 --samFlagExclude 4 --minMappingQuality 3
  samtools index -@ {threads} {output.bam}
  '''

rule CNR_bw:
  input:
    bam = 'CutnRun/sieve/{sample}.bam'
  output:
    bw = 'CutnRun/bws/{sample}.bw'
  threads: 10
  shell:'''
  bamCoverage -b {input.bam} -o {output.bw} -bs 10 -p {threads} --normalizeUsing RPKM
  '''

rule CNR_bw_compare:
  input:
    bam = 'CutnRun/sieve/{chip}.bam',
    allbams = expand('CutnRun/sieve/{sample}.bam', sample=SAMPLES)
  output:
    bw = 'CutnRun/bws/{chip}.log2.bw'
  params:
    ctrl = lambda wildcards: 'CutnRun/sieve/' + chipDict[wildcards.chip] + '.bam'
  threads: 10
  shell:'''
  bamCompare -b1 {input.bam} -b2 {params.ctrl} \
    -p {threads} -o {output.bw} -bs 10
  '''