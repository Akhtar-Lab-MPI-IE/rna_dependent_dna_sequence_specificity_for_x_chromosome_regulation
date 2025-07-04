
with open(snakemake.input.mot,'r') as f:
    with open(snakemake.output.mot, 'w') as o:
        ret = False
        for line in f:
            if line.startswith('MOTIF'):
                name = '>' + line.strip().split(' ')[1]
                o.write(name + '\n')
                ret = True
            elif line.strip() and ret and not line.startswith('letter'):
                o.write(line.strip().replace('  ', '\t') + '\n')
