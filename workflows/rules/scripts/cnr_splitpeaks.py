import pandas as pd

def parseuropa(tsv):
    df = pd.read_table(tsv)
    df['peak'] = df['peak_chr'] + ':' + df['peak_start'].astype(str) + '-' + df['peak_end'].astype(str)
    df = df[
        ['peak', 'peak_chr', 'peak_start', 'peak_end', 'relative_location', 'distance', 'gene_id', 'gene_name', 'feature', 'transcript_biotype']
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
            if r['relative_location'] == 'OverlapStart' and r['feature'] == 'transcript':
                annFeature.append('promoter')
            elif r['relative_location'] == 'OverlapStart' and r['feature'] == 'five_prime_utr':
                annFeature.append('promoter')
            else:
                annFeature.append(r['feature'])
        else:
            annFeature.append('intergenic')
    df['annFeature'] = annFeature
    return(df)


def splitpeaks(annfile, ofx1, ofx2, ofx3, ofnonx1, ofnonx2, ofnonx3):
    _a = parseuropa(annfile)

    # X
    xdf = _a[_a['peak_chr'] == 'X']
    xdf[xdf['annFeature'] == 'intergenic'][['peak_chr', 'peak_start', 'peak_end']].to_csv(
        ofx1, sep='\t', index=None, header=None
    )
    xdf[xdf['annFeature'] == 'promoter'][['peak_chr', 'peak_start', 'peak_end']].to_csv(
        ofx2, sep='\t', index=None, header=None
    )
    xdf[xdf['annFeature'].isin(['exon', 'transcript', 'three_prime_utr', 'five_prime_utr'])][['peak_chr', 'peak_start', 'peak_end']].to_csv(
        ofx3, sep='\t', index=None, header=None
    )
    # Non-X
    nonxdf = _a[_a['peak_chr'] != 'X']
    nonxdf[nonxdf['annFeature'] == 'intergenic'][['peak_chr', 'peak_start', 'peak_end']].to_csv(
        ofnonx1, sep='\t', index=None, header=None
    )
    nonxdf[nonxdf['annFeature'] == 'promoter'][['peak_chr', 'peak_start', 'peak_end']].to_csv(
        ofnonx2, sep='\t', index=None, header=None
    )
    nonxdf[nonxdf['annFeature'].isin(['exon', 'transcript', 'three_prime_utr', 'five_prime_utr'])][['peak_chr', 'peak_start', 'peak_end']].to_csv(
        ofnonx3, sep='\t', index=None, header=None
    )

splitpeaks(
    snakemake.input.peakset,
    snakemake.output.ofx1,
    snakemake.output.ofx2,
    snakemake.output.ofx3,
    snakemake.output.ofnonx1,
    snakemake.output.ofnonx2,
    snakemake.output.ofnonx3
)