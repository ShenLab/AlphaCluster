import pandas as pd
import pickle


def read_fasta(path, ind=0):
    seqs = {}
    seq = ''
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if seq != '':
                    seqs[name] = seq
                    seq = ''
                name = line[1:].strip().split('|')[ind]
            else:
                seq += line.strip()
        if seq != '':
            seqs[name] = seq

        return seqs


def _parse_range(x):
    res = []
    for r in x.split(','):
        if '-' in r[1:]:
            a, b = r.split('-')
        else:
            a, b = r, r
        res.append((a, b))
    return res


def process_pdb_map(r, ensembl_seqs, uniprot_seqs):
    uniprot_ranges = _parse_range(r['uniprot_res'].strip('[]'))
    model_ranges = _parse_range(r['model_res'].strip('[]'))

    check_uniprot_num = 0
    check_model_num = 0
    seq_pdb_set = set()
    '''
    for a, b in uniprot_ranges:
        check_uniprot_num += b - a + 1
        for i in range(a, b + 1):
            seq_pdb_set.add(i)
    for a, b in model_ranges:
        check_model_num += b - a + 1
    assert check_model_num == check_uniprot_num
    '''
    start, end = int(uniprot_ranges[0][0]), int(uniprot_ranges[-1][-1])

    ensembl_seq = ensembl_seqs[r['transcript_id']]
    sub_seq = ensembl_seq[start - 1:end]
    with open(f'./covered/{r["transcript_id"]}.fasta', 'w') as fw:
        fw.write(f'>{r["transcript_id"]}\t{r["uniprot_res"]}\n{sub_seq}\n')

    with open(f'./covered/{r["transcript_id"]}.seq_pdb_map.pickle',
              'wb') as fw:
        pickle.dump(seq_pdb_set, fw)

    return start, end


def main():
    model_df = pd.read_csv('all_models_index.txt', sep='\t')

    uniprot_df = pd.read_csv('./uniprot_map.tsv', sep='\t')
    uniprot_df = uniprot_df.rename(columns={'Entry': 'uniprot_acc'})

    df = pd.merge(model_df, uniprot_df, on='uniprot_acc', how='inner')

    df = df[df['priority'] == 1]

    df['transcript_id'] = df['transcript_id'].apply(lambda x: x.split(',')[0])

    def _filter_canonical(x):
        if type(x) == str and x.split('-')[-1] != '1':
            return False
        return True

    df = df[df['isomap'].apply(_filter_canonical)]

    ensembl_seqs = read_fasta('./all_ensmbl.fasta')
    uniprot_seqs = read_fasta('./all_uniprot.fasta', ind=1)

    def _check_seq(r):
        s1 = ensembl_seqs.get(r['transcript_id'], None)
        s2 = uniprot_seqs.get(r['uniprot_acc'], None)
        if s1 == s2 or len(s1) == len(s2):
            return True

        return False

    df = df[df.apply(_check_seq, axis=1)]
    df.to_csv('uniprot_ensembl_pdb.tsv', index=False, sep='\t')

    info = []
    for idx, r in df.iterrows():
        start, end = process_pdb_map(r, ensembl_seqs, uniprot_seqs)
        info.append((r['transcript_id'], start, end))

    info_df = pd.DataFrame(
        dict(zip(('transcript_id', 'start', 'end'), zip(*info))))

    df2 = pd.merge(df, info_df, on='transcript_id')
    df2.to_csv('./all_ensembl_pdb.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
