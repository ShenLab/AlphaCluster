import time
import json
import argparse
import os
import sys
import logging
import shutil
from datetime import datetime
import glob
import random

import tensorflow as tf
import tensorflow_addons as tfa

from model_attention import ModelAttention
from dataset import build_dataset_2d
from loss import compute_loss

from tensorflow.keras.mixed_precision import experimental as mixed_precision

os.environ['CUDA_VISIBLE_DEVICES'] = '-1'

tf.config.threading.set_intra_op_parallelism_threads(60)
tf.config.threading.set_inter_op_parallelism_threads(60)

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
logging_formatter = logging.Formatter(
    '%(asctime)s - %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
ch = logging.StreamHandler(sys.stdout)
ch.setFormatter(logging_formatter)
logger.addHandler(ch)


def train_single_gpu(config):
    #setup logger
    str_t = datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
    train_dir = f'./res/{str_t}'
    config['train']['train_dir'] = train_dir
    os.makedirs(train_dir)

    fh = logging.FileHandler(f'{train_dir}/train.log')
    fh.setFormatter(logging_formatter)
    logger.addHandler(fh)

    logger.info(json.dumps(config, indent=4))

    batch_size = config['train']['batch_size']
    learning_rate = config['train']['learning_rate']
    input_config = config['input']
    input_base_dir = input_config['base_dir']
    all_files = glob.glob(input_base_dir + '/' + input_config['train'])
    random.shuffle(all_files)
    train_files, validate_files = all_files, all_files

    print('train_file', train_files)

    train_dataset = build_dataset_2d(train_files, batch_size)
    validate_dataset = build_dataset_2d(validate_files, batch_size)

    model = ModelAttention(config['model']['attention'])
    optimizer = tf.keras.optimizers.Adam(learning_rate=learning_rate)

    #metrics
    metric_train_loss = tf.keras.metrics.Mean(name='train_loss')
    metric_test_loss = tf.keras.metrics.Mean(name='test_loss')

    #summary
    train_log_dir = f'{train_dir}/summary/train'
    train_summary_writer = tf.summary.create_file_writer(train_log_dir)

    def _update_histogram_summary():
        with train_summary_writer.as_default():
            for var in model.trainable_variables:
                if 'kernel:' in var.name or 'gamma:' in var.name or 'beta:' in var.name:
                    tf.summary.histogram(var.name,
                                         var,
                                         step=optimizer.iterations)

    def _update_gradient_norm_summary(var, grad):
        with train_summary_writer.as_default():
            for v, g in zip(var, grad):
                if 'kernel:' in v.name or 'gamma:' in v.name or 'beta:' in v.name:
                    tf.summary.scalar(f'gradient_norm/{v.name}',
                                      tf.norm(g, ord='euclidean'),
                                      step=optimizer.iterations)

    @tf.function(input_signature=[validate_dataset.element_spec])
    def test_step(sample):
        ref_aa, alt_aa, feature, label, padding_mask = sample

        logit = model((ref_aa, alt_aa, feature), False, padding_mask)

        loss = compute_loss(label, logit)

        metric_test_loss.update_state(loss)

        pred = tf.math.argmax(logit, axis=-1, output_type=tf.int32)

        mask = tf.math.greater(label, -1.0)

    def test(test_dataset, data_name, epoch):
        metric_test_loss.reset_states()
        for step, sample in enumerate(test_dataset):
            test_step(sample)

        logger.info(
            f'{data_name} aa_num= {metric_test_acc.count.numpy():.0f} acc= {metric_test_acc.result():.3f} loss= {metric_test_loss.result()}'
        )
        return metric_test_loss.result()

    @tf.function(input_signature=[train_dataset.element_spec])
    def train_step(sample):
        ref_aa, alt_aa, feature, label, padding_mask = sample

        with tf.GradientTape() as tape:
            logit = model((ref_aa, alt_aa, feature), True, padding_mask)
            loss = compute_loss(label, logit)

        gradients = tape.gradient(loss, model.trainable_variables)
        optimizer.apply_gradients(zip(gradients, model.trainable_variables))

        metric_train_loss.update_state(loss)
        if optimizer.iterations % 512 == 0:
            _update_gradient_norm_summary(model.trainable_variables, gradients)

        return loss

    EPOCHS = 512
    watch_loss = 10000.0
    watch_epoch = -1
    patience_epochs = 5
    for epoch in range(EPOCHS):
        start = time.time()

        for step, samples in enumerate(train_dataset):
            loss = train_step(samples)
            #model summary
            if optimizer.iterations == 1:
                model.summary(print_fn=logger.info)

            #logging kernel weights
            if (optimizer.iterations + 1) % 512 == 0:
                _update_histogram_summary()

        logger.info(f'Epoch {epoch} Loss {metric_train_loss.result():.4f}')
        metric_train_loss.reset_states()

        model.save(f'{train_dir}/model/epoch-{epoch}')

        #validate and test
        validate_loss = test(validate_dataset, 'validate', epoch)
        if validate_loss < watch_loss:
            watch_loss = validate_loss
            watch_epoch = epoch

        for k, v in config['input']['test'].items():
            test_dataset = build_dataset_2d([v], batch_size)
            test(test_dataset, k, epoch)

        if epoch - watch_epoch == patience_epochs:
            logger.info(f'best_epoch {watch_epoch} min_loss= {watch_loss}')
            break


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str, required=True)
    args = parser.parse_args()

    with open(args.config) as f:
        config = json.load(f)

    train_single_gpu(config)


if __name__ == '__main__':
    main()
