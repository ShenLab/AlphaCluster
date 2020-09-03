import tensorflow as tf

orig_feature_dim = 1 + 1 + 20 + 10 + 192 + 4 + 4
feature_1d_dim = 20 + 10 + 4
feature_3d_dim = 4
used_feature_dim = feature_1d_dim + feature_3d_dim


def build_dataset_2d(input_tfrecord_files, batch_size):
    drop_remainder = False

    feature_description = {
        'label': tf.io.FixedLenFeature([], tf.int64),
        'ref_aa': tf.io.FixedLenFeature([], tf.int64),
        'alt_aa': tf.io.FixedLenFeature([], tf.int64),
        'feature': tf.io.FixedLenFeature([], tf.string),
        'var_id': tf.io.FixedLenFeature([], tf.string),
    }

    def _parser(example_proto):
        parsed = tf.io.parse_single_example(example_proto, feature_description)
        label, ref_aa, alt_aa = parsed['label'], parsed['ref_aa'], parsed[
            'alt_aa']
        var_id = parsed['var_id']

        ref_aa, alt_aa, label = tf.cast(ref_aa, tf.int32), tf.cast(
            alt_aa, tf.int32), tf.cast(label, tf.float32)

        feature = tf.io.decode_raw(parsed['feature'], tf.float32)

        feature = tf.reshape(feature, (-1, orig_feature_dim))

        feature = tf.concat([
            feature[:, 2:32],
            feature[:, 32 + 192:],
        ],
                            axis=-1)
        #feature_1d:20+10+4
        #feature_3d:4

        L = tf.shape(feature)[0]
        padding_mask = tf.zeros((1, L), dtype=tf.float32)

        return ref_aa, alt_aa, feature, label, padding_mask

    dataset = tf.data.TFRecordDataset(input_tfrecord_files)

    options = tf.data.Options()
    options.experimental_threading.max_intra_op_parallelism = 1
    dataset = dataset.with_options(options)

    dataset = dataset.map(_parser, num_parallel_calls=4)

    dataset = dataset.padded_batch(
        batch_size,
        padded_shapes=((), (), (None, used_feature_dim), (), (1, None)),
        padding_values=(tf.cast(-1,
                                tf.int32), tf.cast(-1,
                                                   tf.int32), 0.0, -1.0, 1.0))
    dataset = dataset.prefetch(16)

    return dataset
