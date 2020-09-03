import pandas as pd
import argparse
import pickle
with open('./covered_index.pickle', 'rb') as fr:
    covered_index = pickle.load(fr)

with open('./model_index.pickle', 'rb') as fr:
    model_index = pickle.load(fr)


def build(path):
    output = path.split('.csv')[0] + '.covered.csv'
    df = pd.read_csv(path, sep='\t')

    #df = df[df['transcript_id'].isin(covered_index)]
    header = True
    total = 0
    for tr, g in df.groupby('transcript_id'):
        index = covered_index.get(tr, None)
        if index is None:
            continue
        g = g[g['aa_pos'].isin(index)]
        g['model_start_idx'] = min(index)
        g['model_name'] = model_index.get(tr, None)
        if header:
            g.to_csv(output, index=False, sep='\t')
        else:
            g.to_csv(output, index=False, mode='a', header=False, sep='\t')
        if header:
            header = False
        total += g.shape[0]

    print(df.shape[0], total)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, required=True)
    args = parser.parse_args()
    build(args.input)


if __name__ == '__main__':
    main()
