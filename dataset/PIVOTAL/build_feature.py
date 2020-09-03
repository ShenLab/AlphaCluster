import random
import pandas as pd
import argparse
import util
import traceback
import numpy as np
import json
import os
from multiprocessing import Process
import pickle

random.seed(2020)

compara_dir = '/data/hz2529/zion/missense/compara'
feature_dir = '/data/hz2529/zion/missense/GRCh38/feature'
gnomad_dir = '/data/hz2529/zion/missense/GRCh38_gnomad'
topmed_dir = '/data/hz2529/zion/missense/GRCh38_topmed'
mutation_rate_dir = '/data/hz2529/zion/missense/gnomad-public/split'
compara_dir = '/data/hz2529/zion/missense/GRCh38_gnomad'
output_dir = './feature'


def read_compara(path):
    res = []
    with open(path) as f:
        for line in f:
            seq = line.strip().split('\t')[1]
            seq = np.array([util.aa_index(a) for a in seq])
            res.append(seq)
    return np.array(res).transpose()


def calc_angle(x):
    x = x / 180.0 * np.pi
    return np.stack([np.sin(x), np.cos(x)], axis=-1)


def read_netsurfp2(netsurfp2_path):
    res = []
    with open(netsurfp2_path, 'r') as fr:
        data = json.load(fr)
    names = ['q3_prob', 'rsa', 'interface', 'disorder', 'phi', 'psi']
    for n in names:
        fea = np.array(data[n])
        if n in ['phi', 'psi']:
            fea = calc_angle(fea)
        if fea.ndim == 1:
            fea = np.expand_dims(fea, axis=1)
        res.append(fea)
    res = np.concatenate(res, axis=-1)
    return res


def read_region(region_path):
    res = np.load(region_path)
    return np.nan_to_num(res)


def build_one_transcript(transcript_id, start, end, output_path):
    hhblits_path = f'{feature_dir}/{transcript_id}.hhm.npy'
    netsurfp2_path = f'{feature_dir}/{transcript_id}.netsurfp2.json'
    compara_path = f'{compara_dir}/{transcript_id}.compara'
    region_path = f'{mutation_rate_dir}/{transcript_id}.aa.obs_exp.npy'

    if not os.path.exists(compara_path) or not os.path.exists(
            region_path) or not os.path.exists(netsurfp2_path):
        return

    compara = read_compara(compara_path)
    hhblits = np.load(hhblits_path)[:, 0:21]
    struc = read_netsurfp2(netsurfp2_path)
    region = read_region(region_path)

    compara = compara[start - 1:end]
    hhblits = hhblits[start - 1:end]
    struc = struc[start - 1:end]
    region = region[start - 1:end]
    protein = np.concatenate([hhblits, struc, region, compara], axis=1)
    with open(output_path, 'wb') as fw:
        pickle.dump(protein, fw)


def build(input_path):
    df = pd.read_csv(input_path, sep='\t')
    for idx, r in df.iterrows():
        try:
            build_one_transcript(r['transcript_id'], r['start'], r['end'],
                                 f'{output_dir}/{r["transcript_id"]}.pickle')
        except:
            print('error', r['transcript_id'])
            print(traceback.format_exc())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, required=True)
    args = parser.parse_args()

    build(args.input)
