import os
import glob

def collate_homer(hdir, cutoff):
    '''
    Collates known and unknown motifs from homer, and generates 1 single motif file (homer format)
    directory is supposed to be structured as such, for homer:

    hdir / {some_dir} / homer / knownResults / {knownmotifs}.motif
    hdir / {some_dir} / homer / homerMotifs.all.motifs
    '''
    ph = []
    # Knowns
    _knowns = glob.glob(os.path.join(hdir, '*', 'homer', 'knownResults', '*.motif'))
    for _k in _knowns:
        with open(_k, 'r') as f:
            for line in f:
                ph.append(line.strip())
    # de novo
    _denovo = glob.glob(os.path.join(hdir, '*', 'homer', 'homerMotifs.all.motifs'))
    for _d in _denovo:
        with open(_d, 'r') as f:
            retstatus = False
            for line in f:
                if line.startswith('>'):
                    pval = float(line.strip().split('\t')[5].split(',')[-1].split(':')[1])
                    if pval < cutoff:
                        retstatus = True
                    else:
                        retstatus = False
                if retstatus:
                    ph.append(line.strip())
    return ph

def collate_meme(mdir, mot, cutoff):
    '''
    Collates results from meme-chip.
    Assumes mdir is such that:
    mdir / {some_dir} / meme / combined.meme
    exists for all some_dirs to be collated.
    '''
    _memes = list(glob.glob(os.path.join(mdir, '*', 'meme', 'combined.meme')))
    ph = []
    for _m in _memes:
        # Read significance values for this motif
        sigmotifs = []
        _sumfile = _m.replace('combined.meme', 'summary.tsv')
        if _sumfile.endswith('.tsv'):
            with open(_sumfile, 'r') as sigr:
                for line in sigr:
                    if line.startswith('MOTIF'):
                        continue
                    if not line.strip():
                        continue
                    if line.startswith("#"):
                        continue
                    print(line.strip())
                    if float(line.strip().split('\t')[7]) < cutoff:
                        sigmotifs.append(line.strip().split('\t')[2])

        with open(_m, 'r') as f:
            ret = False
            for line in f:
                if line.startswith('MOTIF'):
                    if line.strip().split(' ')[2].split('-')[0] in sigmotifs:
                        ret = True
                    else:
                        ret = False
                if ret:
                    ph.append(line.strip())
    # Include the extra motif
    with open(mot, 'r') as f:
        ret = False
        for line in f:
            if line.startswith('MOTIF'):
                # No significance cut-off needed here.
                ret = True 
            if ret:
                ph.append(line.strip())
    return ph

def combine(homer, meme, of):
    memeformat = [
        'MEME version 5.5.5',
        '',
        'ALPHABET= ACGT',
        '',
        'strands: + -',
        ''
    ]
    motcount = 0
    for _h in homer:
        if _h.startswith('>'):
            motcount += 1
            if memeformat:
                # Empty line after each motif
                memeformat.append('')
            _ = _h.split('\t')[1]
            motname = f"MOTIF {motcount} {_}"
            memeformat.append(motname)
            memeformat.append('')
            memeformat.append('letter-probability matrix:')
        else:
            memeformat.append(_h.strip())
    for _h in meme:
        if _h.startswith('MOTIF'):
            motcount += 1
            motname = f"MOTIF {motcount} {_h.split(' ')[2]}"
            memeformat.append('')
            memeformat.append(motname)
        else:
            memeformat.append(_h.replace('  ', ''))
    with open(of, 'w') as f:
        for i in memeformat:
            f.write(i + '\n')

combine(
    collate_homer(snakemake.params.hdir, snakemake.params.cutoff),
    collate_meme(snakemake.params.hdir, snakemake.params.extra_mot, snakemake.params.cutoff),
    snakemake.output.mot
)