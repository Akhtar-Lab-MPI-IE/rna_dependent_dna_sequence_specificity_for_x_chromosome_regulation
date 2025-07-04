def ret_regionlabel(_l):
  return ' '.join([i.replace('CNR_MSL1_', '') for i in _l])

rule CNR_computematrix_MSL1:
  input:
    bedfiles = expand('CutnRun/motifs/{peakcat}/peaks.bed', peakcat=peakCategories),
    bw_rep1 = 'CutnRun/bws/CNR_MSL1_rep1.log2.bw',
    bw_rep2 = 'CutnRun/bws/CNR_MSL1_rep2.log2.bw'
  output:
    'CutnRun/deeptools/MSL1peaks.npz'
  threads: 20
  shell:'''
  computeMatrix reference-point -p {threads} \
    -R {input.bedfiles} \
    -S {input.bw_rep1} {input.bw_rep2} \
    -o {output} \
    -b 5000 -a 5000 \
    --missingDataAsZero \
    --referencePoint center
  '''

rule CNR_plotheatmap_MSL1:
  input:
    'CutnRun/deeptools/MSL1peaks.npz'
  output:
    'figures/cutnrun_heatmap_MSL1.pdf'
  params:
    rlabels = ret_regionlabel(peakCategories)
  threads: 2
  shell:'''
  plotHeatmap -m {input} -o {output} \
    --regionsLabel {params.rlabels} \
    --refPointLabel center \
    --samplesLabel MSL1_rep1 MSL1_rep2 \
    --whatToShow "heatmap and colorbar"
  '''
