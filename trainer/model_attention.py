import tensorflow as tf
import tensorflow_addons as tfa
from model import ModelBase

import numpy as np


def get_angles(pos, i, d_model):
    angle_rates = 1 / np.power(10000, (2 * (i // 2)) / np.float32(d_model))
    return pos * angle_rates


def positional_encoding(position, d_model):
    angle_rads = get_angles(
        np.arange(position)[:, np.newaxis],
        np.arange(d_model)[np.newaxis, :], d_model)

    angle_rads[:, 0::2] = np.sin(angle_rads[:, 0::2])
    angle_rads[:, 1::2] = np.cos(angle_rads[:, 1::2])

    pos_encoding = angle_rads[np.newaxis, ...]

    return tf.cast(pos_encoding, dtype=tf.float32)


def positional_encoding_2d(position, d_model):
    x = np.arange(position)[:, np.newaxis]
    x = np.tile(x, (1, position))
    y = x.transpose()
    z = np.abs(x - y)
    angle_rads = get_angles(z[:, :, np.newaxis],
                            np.arange(d_model)[np.newaxis, :], d_model)
    # apply sin to even indices in the array; 2i
    angle_rads[:, :, 0::2] = np.sin(angle_rads[:, :, 0::2])

    # apply cos to odd indices in the array; 2i+1
    angle_rads[:, :, 1::2] = np.cos(angle_rads[:, :, 1::2])

    pos_encoding = angle_rads[np.newaxis, ...]

    return tf.cast(pos_encoding, dtype=tf.float32)


class MultiHeadDot(tf.keras.layers.Layer):
    def __init__(self, d_model, num_heads):
        super(MultiHeadDot, self).__init__()

        self.kernel = self.add_weight(name='dot_kernel',
                                      shape=(d_model, ),
                                      dtype=tf.float32,
                                      initializer='glorot_uniform',
                                      trainable=True)
        self.num_heads = num_heads
        self.d_model = d_model

    def call(self, inputs):

        x = tf.math.multiply(inputs, self.kernel)
        shape = tf.shape(x)
        x = tf.reshape(x, (shape[0], shape[1], self.num_heads, -1))

        x = tf.reduce_sum(x, axis=-1)  #(batch_size, neighbour_size, num_heads)
        return x

    def get_config(self):
        config = super(MultiHeadDot, self).get_config()
        config.update({"d_model": self.d_model, "num_heads": self.num_heads})
        return config


class MultiHeadAttention(tf.keras.layers.Layer):
    def __init__(self, d_model, num_heads):
        super(MultiHeadAttention, self).__init__()
        self.num_heads = num_heads
        self.d_model = d_model

        assert d_model % self.num_heads == 0

        self.depth = d_model // num_heads

        self.wq = tf.keras.layers.Dense(d_model)
        self.wk = tf.keras.layers.Dense(d_model)
        self.wv = tf.keras.layers.Dense(d_model)
        self.w2d = tf.keras.layers.Dense(d_model)

        self.mdot = MultiHeadDot(d_model, num_heads)

        self.attention_activation = tf.keras.activations.tanh

        #self.linear_proj = tf.keras.layers.Dense(d_model)

    def call(self, x, mask=None):
        x1d, x3d = x[:, :34], x[:, 34:]
        shape = tf.shape(x1d)
        batch_size, neighbour_size = shape[0], shape[1]

        q, k, v = x1d[:, :1], x1d, x1d

        #set query only for center position
        q = self.wq(q)  # (batch_size, 1, d_model)

        k = self.wk(k)  # (batch_size, neighbour_size, d_model)
        v = self.wv(v)  # (batch_size, neighbour_size, d_model)
        x3d = self.w3d(x3d)  #(batch_size, neighbour_size, d_model)

        #additive attention
        attention_logits = self.mdot(self.attention_activation(
            q + k + x3d))  #(batch_size, neighbour_size, num_heads)

        #split heads for attention logits
        attention_logits = tf.transpose(
            attention_logits,
            perm=[0, 2, 1])  #(batch_size, num_heads, neighbour_size)

        scaled_attention_logits = attention_logits / tf.math.sqrt(
            tf.cast(self.depth, tf.float32))

        #mask, (batch_size, 1, neighbour_size)
        if mask is not None:
            scaled_attention_logits += (mask * -1e9)

        attention_weights = tf.nn.softmax(
            scaled_attention_logits,
            axis=-1)  # (batch_size, num_heads,  neighbour_size)

        #split heads for values
        v = tf.reshape(
            v, (batch_size, -1, num_heads,
                depth))  #(batch_size, neighbour_size, num_heads, depth)
        v = tf.transpose(
            x, perm=[0, 2, 1,
                     3])  # (batch_size, num_heads, neighbour_size, depth)

        scaled_attention = tf.squeeze(
            tf.matmul(attention_weights[:, :, tf.newaxis, :],
                      v), axis=2)  # (batch_size, num_heads,  depth)

        concat_attention = tf.reshape(
            scaled_attention,
            (batch_size, self.d_model))  # (batch_size, d_model)

        return output


class ModelAttention(ModelBase):
    def __init__(self, model_config):
        super(ModelAttention, self).__init__()

        num_layers = model_config['num_layers']
        d_model = model_config['d_model']
        num_heads = model_config['num_heads']
        dff = model_config['dff']
        rate = model_config['rate']
        maximum_position_encoding = model_config['maximum_position_encoding']

        self.pos_encoding_2d = positional_encoding_2d(
            maximum_position_encoding, d_model)

        self.mha = MultiHeadAttention(d_model, num_heads)

        self.proj_1d = tf.keras.layers.Dense(d_model, activation='relu')
        self.proj_3d = tf.keras.layers.Dense(d_model, activation='relu')

        layer_w = [d_model, d_model // 2]
        self.dense_layers = [
            tf.keras.layers.Dense(w, activation='relu') for w in layer_w
        ]

        self.logit_layer = tf.keras.layers.Dense(3)

    def call(self, inputs, training=False, mask=None):
        ref_aa, alt_aa, feature = inputs
        feature_1d, feature_3d = feature[:, :, :34], feature[:, :, 34:]

        feature_1d = self.proj_1d(feature_1d)  #(batch_size, len, model_dim)
        feature_3d = self.proj_2d(feature_3d)  #(batch_size, len, model_dim)

        context_3d = self.mha(feature_1d, feature_3d, mask)

        ref_aa = tf.one_hot(tf.cast(ref_aa, tf.int32),
                            depth=20,
                            dtype=tf.float32)
        alt_aa = tf.one_hot(tf.cast(alt_aa, tf.int32),
                            depth=20,
                            dtype=tf.float32)

        x = tf.concat([ref_aa, alt_aa, context_3d], axis=-1)

        for n in self.dense_layers:
            x = n(x)

        x = self.logit_layer(x)

        return x
