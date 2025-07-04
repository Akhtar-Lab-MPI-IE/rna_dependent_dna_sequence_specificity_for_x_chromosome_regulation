SAMPLES, = glob_wildcards("fqs/FLASH/{sample}_R1.fastq.gz")

CTRLDIC = {
  'Dme_Clone8_MSL1_rep1': 'Dme_Clone8_IgG_rep1',
  'Dme_Clone8_MSL1_rep2': 'Dme_Clone8_IgG_rep2',
  'Dme_Clone8_MSL1_rep3': 'Dme_Clone8_IgG_rep3',
  'Dme_S2_MSL1_rep1': 'Dme_S2_IgG_rep1',
  'Dme_S2_MSL1_rep2': 'Dme_S2_IgG_rep2',
  'Dme_S2_MSL1_rep3': 'Dme_S2_IgG_rep3',
  'Dme_MLE_male_rep1': 'Dme_WOR_male',
  'Dme_MLE_male_rep2': 'Dme_WOR_male',
  'Dme_MLE_female_rep1': 'Dme_WOR_female',
  'Dme_MLE_female_rep2': 'Dme_WOR_female',
  'Dme_MSL2_male_rep1': 'Dme_WOR_male',
  'Dme_MSL2_male_rep2': 'Dme_WOR_male',
  'Dme_MSL2_male_rep3': 'Dme_WOR_male',
  'Dme_MSL2_female_rep1': 'Dme_WOR_female',
  'Dme_MSL2_female_rep2': 'Dme_WOR_female',
  'Dme_MSL2_female_rep3': 'Dme_WOR_female'
}

MALE_FEMALEDIC = {
    'Dme_MSL2_male_rep1': 'Dme_MSL2_female_rep1',
    'Dme_MSL2_male_rep2': 'Dme_MSL2_female_rep2',
    'Dme_MSL2_male_rep3': 'Dme_MSL2_female_rep3',
    'Dme_MLE_male_rep1': 'Dme_MLE_female_rep1',
    'Dme_MLE_male_rep2': 'Dme_MLE_female_rep2'
}

CORDIC = {
  'MSL1_C8': ['Dme_Clone8_MSL1_rep1', 'Dme_Clone8_MSL1_rep2', 'Dme_Clone8_MSL1_rep3', 'Dme_Clone8_IgG_rep1', 'Dme_Clone8_IgG_rep2', 'Dme_Clone8_IgG_rep3'],
  'MSL1_S2': ['Dme_S2_MSL1_rep1', 'Dme_S2_MSL1_rep2', 'Dme_S2_MSL1_rep3', 'Dme_S2_IgG_rep1', 'Dme_S2_IgG_rep2', 'Dme_S2_IgG_rep3'],
  'MLE': ['Dme_MLE_male_rep1', 'Dme_MLE_male_rep2', 'Dme_MLE_female_rep1', 'Dme_MLE_female_rep2'],
  'MSL2': ['Dme_MSL2_male_rep1' ,'Dme_MSL2_male_rep2', 'Dme_MSL2_male_rep3', 'Dme_MSL2_female_rep1', 'Dme_MSL2_female_rep2', 'Dme_MSL2_female_rep3']
}

include: 'rules/flash_align.smk'
include: 'rules/flash_dewseq.smk'

rule all:
  input:
    'FLASH/multiqc/multiqc_report.html',
    expand('FLASH/peakannotation/{ipsample}_finalhits.txt', ipsample = CTRLDIC.keys()),
    'figures/FLASH_annotation.pdf',
    expand('FLASH/bw/{ipsample}.log2.bw', ipsample = CTRLDIC.keys()),
    expand('FLASH/plotcorr/{CORR}_plotcorrRPKM.pdf', CORR=CORDIC.keys()),
    # bw male vs female
    expand('FLASH/bw_mf/{malef}.mvsf.log2.bw', malef = MALE_FEMALEDIC.keys()),
    # dewseq
    'FLASH/dewseq/FLASH_enrichment.pdf'
