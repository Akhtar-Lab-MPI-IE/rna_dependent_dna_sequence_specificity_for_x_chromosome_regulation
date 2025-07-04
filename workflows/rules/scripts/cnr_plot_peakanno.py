import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

cdic = {
    "three_prime_utr":"#ff7f0e",
    "five_prime_utr":"#9467bd",
    "first_exon":"#2ca02c",
    "intergenic":"#7f7f7f",
    "promoter":"#1f77b4",
    "first_intron":"#d62728",
    "other_exon": "#8c564b",
    "other_intron": "#e377c2"
}

def parsegtf(gtf):
    '''
    create genedic
    genedic[geneid] = {
        txn = (start, end)
        exons = [(start, end), (start, end), ...]
    }
    '''
    genedic = {}
    with open(gtf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            geneid = None
            l = line.strip().split('\t')
            for i in l[8].split(';'):
                if 'gene_id' in i:
                    geneid = i.split(' ')[-1].replace('"', '').strip()
            if geneid:
                if geneid not in genedic:
                    genedic[geneid] = {
                        'txn': None,
                        'exons': []
                    }
                if l[2] == 'transcript':
                    start = int(l[3])
                    end = int(l[4])
                    if not genedic[geneid]['txn']:
                        genedic[geneid]['txn'] = (start, end)
                    else:
                        if start < genedic[geneid]['txn'][0]:
                            genedic[geneid]['txn'] = (start, genedic[geneid]['txn'][1])
                        if end > genedic[geneid]['txn'][1]:
                            genedic[geneid]['txn'] = (genedic[geneid]['txn'][0], end)
                elif l[2] == 'exon':
                    start = int(l[3])
                    end = int(l[4])
                    genedic[geneid]['exons'].append((start, end))
    for i in genedic:
        genedic[i]['exons'] = sorted(genedic[i]['exons'], key=lambda x: x[0])
        # define first intron
        if len(genedic[i]['exons']) > 1:
            genedic[i]['firstintron'] = (genedic[i]['exons'][0][1], genedic[i]['exons'][1][0])
    return genedic

def parseuropa(tsv, genedic):
    df = pd.read_table(tsv)
    df['peak'] = df['peak_chr'] + ':' + df['peak_start'].astype(str) + '-' + df['peak_end'].astype(str)
    df = df[
        ['peak', 'peak_chr', 'peak_start', 'peak_end', 'relative_location', 'gene_id', 'gene_name', 'feature', 'transcript_biotype']
    ]
    # Parse uropa df to combine feature and rel_loc
    annFeature = []
    for i,r in df.iterrows():
        if r['relative_location'] in [
            'OverlapStart',
            'OverlapEnd',
            'PeakInsideFeature',
            'FeatureInsidePeak'
        ]:
            annofeat = None
            if r['relative_location'] == 'OverlapStart' and r['feature'] == 'transcript':
                annFeature.append('promoter')
            elif r['relative_location'] == 'OverlapStart' and r['feature'] == 'five_prime_utr':
                annFeature.append('promoter')
            else:
                # Differentiate between exons and introns, but only if feature is transcript or exon.
                if r['feature'] in ['transcript', 'exon']:
                    geneid = r['gene_id']
                    assert geneid in genedic, f"Gene {geneid} not found in genedic.."
                    start = r['peak_start']
                    end = r['peak_end']
                    if r['feature'] == 'transcript':
                        # Since feature is transcript, it should not be in an exon.
                        for exon in genedic[geneid]['exons']:
                            assert start <= exon[1] or start >= exon[1], f"Feature transcript, but found in exon: {r} <-> {exon}.."
                        assert 'firstintron' in genedic[geneid], f"No first intron annotated, but feature is transcript.. {r} <-> {genedic[geneid]}"
                        if start >= genedic[geneid]['firstintron'][0] and start <= genedic[geneid]['firstintron'][1]:
                            annofeat = 'first_intron'
                        elif end >= genedic[geneid]['firstintron'][0] and end <= genedic[geneid]['firstintron'][1]:
                            annofeat = 'first_intron'
                        else:
                            annofeat = 'other_intron'
                    elif r['feature'] == 'exon':
                        if start >= genedic[geneid]['exons'][0][0] and start <= genedic[geneid]['exons'][0][1]:
                            annofeat = 'first_exon'
                        elif end >= genedic[geneid]['exons'][0][0] and end <= genedic[geneid]['exons'][0][1]:
                            annofeat = 'first_exon'
                        else:
                            annofeat = 'other_exon'
                    annFeature.append(annofeat)
                else:
                    annFeature.append(r['feature'])
        else:
            annFeature.append('intergenic')
    df['annFeature'] = annFeature
    return(df)

def plotpeaks(MSL1, genes, ofile):
    # Parse gtfs for exon - intron definitions
    genedic = parsegtf(genes)
    

    MSL1 = parseuropa(MSL1, genedic)
    
    fig, ax = plt.subplots(ncols=1, nrows=3, figsize=(6,12))

    ######################################### Row1 ######################################### 
    ax[0].pie(
        [len(MSL1[MSL1['peak_chr'] == 'X']), len(MSL1[MSL1['peak_chr'] != 'X'])],
        labels=['X peaks', 'non-X peaks'],
        explode=[0.1,0],
        colors=['#1f77b4', '#7f7f7f'],
        autopct=lambda p: f"{p:.1f}%"
    )
    ax[0].set_title('MSL1')

    ######################################### Row2 ######################################### 

    sns.barplot(
        hue=(MSL1[MSL1['peak_chr'] == 'X']['annFeature'].value_counts()/len(MSL1[MSL1['peak_chr'] == 'X'])).index,
        y=(MSL1[MSL1['peak_chr'] == 'X']['annFeature'].value_counts()/len(MSL1[MSL1['peak_chr'] == 'X'])),
        palette=cdic,
        ax=ax[1]
    )
    ax[1].set_title(f"X peaks: {len(MSL1[MSL1['peak_chr'] == 'X'])} peaks")
    ax[1].set_xlabel("")
    ax[1].set_ylabel("Fraction")
    ax[1].tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False) 
    handles, labels = ax[1].get_legend_handles_labels()
    ax[1].legend(handles=handles, labels=labels)

    ######################################### Row3 ######################################### 

    sns.barplot(
        hue=(MSL1[MSL1['peak_chr'] != 'X']['annFeature'].value_counts()/len(MSL1[MSL1['peak_chr'] != 'X'])).index,
        y=(MSL1[MSL1['peak_chr'] != 'X']['annFeature'].value_counts()/len(MSL1[MSL1['peak_chr'] != 'X'])),
        palette=cdic,
        ax=ax[2],
        legend=False
    )
    ax[2].set_title(f"non-X peaks: {len(MSL1[MSL1['peak_chr'] != 'X'])} peaks")
    ax[2].set_xlabel("")
    ax[2].set_ylabel("Fraction")
    ax[2].tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False) 

    fig.savefig(ofile, dpi=300)

plotpeaks(
    snakemake.input.MSL1,
    snakemake.params.gene,
    snakemake.output[0]
)