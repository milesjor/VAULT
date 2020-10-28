import re
import io
import pandas as pd


def vcf_to_matrix():
    vcf = "/Users/bic/Desktop/mtseq_paper_data/mutation_95p/nsg_3.vcf"
    out_file = "/Users/bic/Desktop/mtseq_paper_data/mutation_95p/heatmap_allvaf/nsg_3.snv.pos.matrix.txt"

    def read_vcf(path):
        with open(path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]
        return pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                   'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t'
        ).rename(columns={'#CHROM': 'CHROM'})

    pd_vcf = read_vcf(vcf)

    # # below is for troubleshooting the bad line
    # print(pd_vcf.iloc[[31330], 0:4])
    # for row in range(31000, len(pd_vcf)):
    #     print("row is ", row)
    #     # print(pd_vcf.at[row, 'POS'])
    #     # print(pd_vcf.iloc[[row], 0:4])
    #     pd_vcf.iloc[[row], [1]] = pd_vcf.iloc[[row], [1]].astype(str).astype(int)

    pd_vcf["POS"] = pd.to_numeric(pd_vcf["POS"])

    count_pos = pd_vcf.POS.value_counts().reset_index().rename(columns={'POS': 'count', 'index': 'POS'})
    pos_sort = count_pos.sort_values("POS").reset_index(drop=True)

    samples = pd_vcf.ID.value_counts().reset_index().rename(columns={'ID': 'count', 'index': 'name'})
    # samples_sort = samples.POS.sort_values().reset_index(drop=True)

    new_df = pd.DataFrame(0, columns=pos_sort.POS, index=samples.name)

    print(new_df)
    with open(vcf, 'r') as infile:
        for line in infile:
            if not line.startswith("#"):
                line = re.split(r'\t', line)
                name = line[2]
                pos = int(line[1])
                new_df.at[name, pos] = 1
    print(new_df)
    
    new_df.to_csv(out_file, header=True, index=True, sep='\t', mode='w')


vcf_to_matrix()
