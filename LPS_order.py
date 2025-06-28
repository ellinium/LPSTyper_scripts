
#E. coli LPS gene order, starts with waa
#Use after converting GBK file to CSV (genbank2CSV_full.py)
import pandas as pd


def main():
    file_name = 'NCTC122.csv'

    df = pd.read_csv(file_name)
    df = df.fillna('missing')

    gene_df = df.loc[df['GeneID'].str.contains('waa')]
    gene_df = gene_df.sort_values(by=['CDS_start'])

    print(str(len(gene_df)))
    lps_seq = ''
    for gene  in gene_df['GeneID']:
        lps_seq = lps_seq + gene + '->'
    print(lps_seq)

if __name__ == '__main__':
    main()