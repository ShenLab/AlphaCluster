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
    s = set()
    for a, b in res:
        for i in range(int(a), int(b) + 1):
            s.add(i)
    return s


def main():
    df = pd.read_csv('../PIVOTAL/all_ensembl_pdb.tsv', sep='\t')

    covered_index = {}
    model_dict = {}
    for idx, r in df.iterrows():
        covered_index[r['transcript_id']] = _parse_range(
            r['uniprot_res'].strip('[]'))
        model_dict[r['transcript_id']] = r['pdb_path'].split('.')[0]

    with open('covered_index.pickle', 'wb') as fw:
        pickle.dump(covered_index, fw)

    with open('model_index.pickle', 'wb') as fw:
        pickle.dump(model_dict, fw)


if __name__ == '__main__':
    main()
