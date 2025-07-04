rule FLASH_dewflatten:
  output:
    'FLASH/dewseq/annotation/genes.bed'
  params:
    gff = config['gff']
  conda: 'envs/dewseq.yml'
  shell:'''
  htseq-clip annotation --gff {params.gff} -o {output} \
    -u gene_id -n Name -t biotype --splitExons --unsorted
  '''

rule FLASH_dewslide:
  input:
    flat = 'FLASH/dewseq/annotation/genes.bed',
  output:
    'FLASH/dewseq/annotation/genes.slidingwindows.txt'
  conda: 'envs/dewseq.yml'
  shell:'''
  htseq-clip createSlidingWindows -i {input.flat} -o {output} -w 150 -s 80
  '''

rule FLASH_dewmappingfile:
  input:
    'FLASH/dewseq/annotation/genes.slidingwindows.txt'
  output:
    'FLASH/dewseq/annotation/genes.map2id.txt'
  conda: 'envs/dewseq.yml'
  shell:'''
  htseq-clip mapToId -a {input} -o {output}
  '''

rule FLASH_sites:
  input:
    'FLASH/dedup/{sample}.bam'
  output:
    'FLASH/dewseq/sites/{sample}.sites.txt'
  threads: 10
  conda: 'envs/dewseq.yml'
  shell:'''
  htseq-clip extract -i {input} -o {output} --cores {threads} \
    --mate 1 --site s
  '''

rule FLASH_count:
  input:
    sites = 'FLASH/dewseq/sites/{sample}.sites.txt',
    windows = 'FLASH/dewseq/annotation/genes.slidingwindows.txt'
  output:
    counts = 'FLASH/dewseq/counts/{sample}.counts.txt'
  conda: 'envs/dewseq.yml'
  shell:'''
  htseq-clip count -i {input.sites} -o {output.counts} \
    -a {input.windows}
  '''

rule FLASH_countmatrix:
  input:
    counts = expand('FLASH/dewseq/counts/{sample}.counts.txt', sample=SAMPLES),
    map2id = 'FLASH/dewseq/annotation/genes.map2id.txt'
  output:
    'FLASH/dewseq/countmatrix.txt'
  params:
    cdir = 'FLASH/dewseq/counts/'
  conda: 'envs/dewseq.yml'
  shell:'''
  htseq-clip createMatrix -i {params.cdir} -o {output} -e counts.txt
  '''

rule FLASH_DEWde:
  input:
    countData = 'FLASH/dewseq/countmatrix.txt',
    annotationData = 'FLASH/dewseq/annotation/genes.map2id.txt'
  output:
    'FLASH/dewseq/DEWde.html'
  params:
    C8out = 'FLASH/dewseq/C8.tsv',
    S2out = 'FLASH/dewseq/S2.tsv',
    MLEout = 'FLASH/dewseq/MLE.tsv',
    MSL2out = 'FLASH/dewseq/MSL2.tsv',
    MLEmalefemaleout = 'FLASH/dewseq/MLE_male_female.tsv',
    MSL2malefemaleout = 'FLASH/dewseq/MSL2_male_female.tsv'
  conda: 'envs/dewseq.yml'
  script:
    'scripts/dewseq.Rmd'

rule FLASH_plotDEW:
  input:
    'FLASH/dewseq/DEWde.html'
  output:
    opng = 'FLASH/dewseq/FLASH_enrichment.pdf'
  params:
    C8 = 'FLASH/dewseq/C8.tsv',
    S2 = 'FLASH/dewseq/S2.tsv',
    MLE = 'FLASH/dewseq/MLE.tsv',
    MSL2 = 'FLASH/dewseq/MSL2.tsv',
    MLEf = 'FLASH/dewseq/MLE_male_female.tsv',
    MSL2f = 'FLASH/dewseq/MSL2_male_female.tsv',
    FDRcut = 0.1
  script:
    'scripts/dewseq_plot_enrich.py'
  
    