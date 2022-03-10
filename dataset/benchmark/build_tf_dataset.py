import random
import pandas as pd
import argparse
import traceback
import numpy as np
import os, sys
from multiprocessing import Process
import pickle

import tensorflow as tf

sys.path.insert(1, '../PIVOTAL')
import util

random.seed(2020)

feature_1d_dir = '../PIVOTAL/feature'
feature_3d_dir = '../PIVOTAL/structure'


def _bytes_feature(value):
    return tf.train.Feature(bytes_list=tf.train.BytesList(value=[value]))


def _int64_feature(value):
    return tf.train.Feature(int64_list=tf.train.Int64List(value=[value]))


def _build_one_var(r, feature_1d, feature_3d):
    L = feature_1d.shape[0]
    aa_pos, ref_aa, alt_aa, label = r['aa_pos'], util.aa_index(
        r['ref_aa']), util.aa_index(r['alt_aa']), r['target']
    var_id, model_start_idx = r['var'], r['model_start_idx']

    feature = []

    protein_feature_idx = list(range(0, 31)) + list(range(35, 35 + 192))

    def _get_one_pos(other_idx, f3d):
        center = other_idx - model_start_idx
        protein_1d = feature_1d[center, protein_feature_idx]
        region_w = 9 // 2
        region_start, region_end = max(0, center - region_w), min(
            L, center + region_w + 1)
        region_1d = np.sum(feature_1d[region_start:region_end, 31:35], axis=0)

        f = np.concatenate([
            np.array([np.abs(other_idx - aa_pos)]), protein_1d, region_1d, f3d
        ],
                           axis=-1)
        feature.append(f)

    _get_one_pos(aa_pos, feature_3d[(str(aa_pos), str(aa_pos))])

    dis_thres = 12.0
    sep_thres = 6.0
    for other_idx in range(L):
        v = feature_3d.get((str(aa_pos), str(other_idx)), None)
        if v is None:
            continue
        if v[0] < dis_thres and np.abs(other_idx - aa_pos) >= sep_thres:
            _get_one_pos(other_idx, v)

    feature = np.array(feature)

    tf_feature = {}
    tf_feature['label'] = _int64_feature(int(label))
    tf_feature['ref_aa'] = _int64_feature(ref_aa)
    tf_feature['alt_aa'] = _int64_feature(alt_aa)
    tf_feature['feature'] = _bytes_feature(
        feature.astype(np.float32).tobytes())
    tf_feature['var_id'] = _bytes_feature(var_id.encode('utf-8'))

    example = tf.train.Example(features=tf.train.Features(feature=tf_feature))

    return example


def build_one_transcript(df, tf_writer):
    transcript_id = df.iloc[0]['transcript_id']
    model_name = df.iloc[0]['model_name']

    feature_1d_path = f'{feature_1d_dir}/{transcript_id}.pickle'
    feature_3d_path = f'{feature_3d_dir}/{model_name}.pickle'

    if not os.path.exists(feature_1d_path) or not os.path.exists(
            feature_3d_path):
        return 0, 0

    with open(feature_1d_path, 'rb') as f1, open(feature_3d_path, 'rb') as f3:
        feature_1d = pickle.load(f1)
        feature_3d = pickle.load(f3)

    pos_num, neg_num = 0, 0

    for idx, r in df.iterrows():
        try:
            example = _build_one_var(r, feature_1d, feature_3d)

            if example is not None:
                tf_writer.write(example.SerializeToString())
            if r['target'] == 0:
                neg_num += 1
            else:
                pos_num += 1
        except:
            traceback.print_exc()

    return pos_num, neg_num


def build_one_thread(curated, output):
    pos_num = 0
    neg_num = 0
    with tf.io.TFRecordWriter(f'{output}.tfrec') as writer:
        for transcript_name, df in curated.groupby('transcript_id'):
            try:
                p, n = build_one_transcript(df, writer)
                pos_num += p
                neg_num += n
            except:
                traceback.print_exc()
    print(f'thread pos_num= {pos_num} neg_num= {neg_num}')


def build(path, output, cpu, af):
    df = pd.read_csv(path, sep='\t')
    df = df[df['gnomad_af'] < af]
    print('target')
    print(df['target'].value_counts())
    print('transcript num')
    print(df['transcript_id'].nunique())

    transcript_list = list(df['transcript_id'].unique())

    random.shuffle(transcript_list)

    if cpu <= 1:
        build_one_thread(df, output)
    else:
        num_each = int((len(transcript_list) - 1) / cpu) + 1
        pool = []
        for idx in range(cpu):
            start = idx * num_each
            end = start + num_each
            if idx == cpu - 1:
                end = len(transcript_list)
            p = Process(target=build_one_thread,
                        args=(df[df['transcript_id'].isin(
                            transcript_list[start:end])], f'{output}_{idx}'))
            pool.append(p)
        for p in pool:
            p.start()
        for p in pool:
            p.join()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, required=True)
    parser.add_argument('--output', type=str, required=True)
    parser.add_argument('--cpu', type=int, default=1)
    parser.add_argument('--af', type=float, default=1.0)
    args = parser.parse_args()

    build(args.input, args.output, args.cpu, args.af)
