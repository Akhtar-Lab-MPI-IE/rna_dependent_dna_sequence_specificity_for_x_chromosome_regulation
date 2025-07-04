from scripts.workflow_functions import calcHash, runSmk, shipResults, lineprinter
from pathlib import Path
import sys
from rich import print

# CONFIG
CONFIGFILE = Path(__file__).parents[0] / 'conf' / 'smk_config.yml'
RESULTSDIR = Path(__file__).parents[0] / 'results'
SMK_PROFILE = 'snakemake_profile'
WDIR = 'WDIR'

# RESULTS
RFS = [
    ('figures/*.pdf', 'figures'),
    ('CutnRun/bws/*bw', 'bws_CutnRun'),
    ('RNA/*html', 'rna'),
    ('RNA/*tsv', 'rna'),
    ('RNA/*pdf', 'rna'),
    ('CutnRun/motifs/*/meme/meme-chip.html', 'motif', 2),
    ('CutnRun/motifs/*/homer/*html', 'motif', 2),
    ('FLASH/bw/*bw', 'bws_FLASH'),
    ('FLASH/bw_mf/*bw', 'bws_FLASH'),
    ('FLASH/dewseq/*tsv', 'dewseq_FLASH'),
    ('FLASH/dewseq/FLASH_enrichment.pdf', 'dewseq_FLASH')
]

# WORKFLOWS
WFS = ['rna.smk', 'cutnrun.smk', 'flash.smk']
rets = []
for _wf in WFS:
    rets.append(
        runSmk(
            Path('workflows', _wf),
            CONFIGFILE,
            WDIR,
            SMK_PROFILE
        )
    )

if rets.count(0) == len(rets):
    tdic = shipResults(RFS, WDIR, RESULTSDIR)
    lineprinter()
    print(f"Results in {RESULTSDIR}:")
    print(tdic)