import random
import pandas as pd
import argparse
import os

random.seed(2020)


def get_curated(repeat, output):
    af_thres = 1e-4
    protein_len_thres = 10000
    protein_len_thres = 5000

    asd = pd.read_csv('./data/GRCh38_ASD.covered.csv', sep='\t')
    chd = pd.read_csv('./data/GRCh38_CHD.covered.csv', sep='\t')
    control = pd.read_csv('./data/GRCh38_Control.covered.csv', sep='\t')

    hgmd = pd.read_csv('./train/HGMD/GRCh38_HGMD.covered.csv', sep='\t')
    clinvar = pd.read_csv('./train/ClinVar/GRCh38_ClinVar.covered.csv',
                          sep='\t')
    uniprot = pd.read_csv('./train/UniProt/GRCh38_UniProt.covered.csv',
                          sep='\t')

    test = pd.concat([asd, chd, control], axis=0)
    test_set = set(list(test['var']))

    df = pd.concat([clinvar, hgmd, uniprot], axis=0)

    pos_id = set(df[df['target'] == 1]['protein_var'])
    neg_id = set(df[df['target'] == 0]['protein_var'])
    over = pos_id & neg_id
    print('number of overlap between positives and negatives', len(over))
    df = df[df['protein_var'].apply(lambda x: x not in over)]
    df = df.drop_duplicates(['var'])
    print('Filter: inconsistent variants.', df.shape[0])

    df = df[df['gnomad_af'] < af_thres]
    print(f'Filter: AF < {af_thres}.', df.shape[0])
    print('Curated variants')
    print(df['target'].value_counts())

    curated = df
    curated = curated[~curated.isin(test_set)]
    curated_set = set(list(curated['var']))
    used_var_set = curated_set | test_set

    random = pd.read_csv('./train/DiscovEHR/GRCh38_DiscovEHR.covered.csv',
                         sep='\t')
    random = random[random['gnomad_af'] < af_thres]
    random['target'] = 0

    random = random[~random['var'].isin(used_var_set)]
    print('random', random.shape[0])

    curated = curated[~curated['genename'].isin(['BRCA1', 'PTEN', 'TP53'])]
    random = random[~random['genename'].isin(['BRCA1', 'PTEN', 'TP53'])]

    curated_pos_num = curated[curated['target'] == 1].shape[0]
    curated_neg_num = curated[curated['target'] == 0].shape[0]
    random_num = curated_pos_num * 2 - curated_neg_num

    for i in range(repeat):
        print('repeat', i, flush=True)
        random_num = min(random_num, random.shape[0])
        neg = random.sample(random_num, random_state=2020 + i)
        res = pd.concat([curated, neg], axis=0)

        print('Filter out the variants in testing genes.', res.shape[0])
        print(res['target'].value_counts())
        res.to_csv(f'{args.output}/r{i}.csv', sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--repeat', type=int, required=True)
    parser.add_argument('--output', type=str, required=True)
    args = parser.parse_args()
    get_curated(args.repeat, args.output)
