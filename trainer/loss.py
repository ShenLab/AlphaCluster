import tensorflow as tf


def compute_per_example_loss(label, logit):
    if tf.shape(label)[0] == 0:
        return tf.cast([], dtype=logit.dtype)

    mask = tf.cast(tf.greater(label, -1.0), dtype=logit.dtype)
    label = tf.one_hot(tf.cast(label, tf.int32), depth=3, dtype=logit.dtype)

    loss = tf.nn.softmax_cross_entropy_with_logits(label, logit)

    loss *= mask

    per_example_loss = tf.reduce_sum(loss, axis=1) / tf.reduce_sum(mask,
                                                                   axis=1)

    return per_example_loss


def compute_loss(label, logit):

    mask = tf.cast(tf.greater(label, -1.0), dtype=logit.dtype)
    label = tf.one_hot(tf.cast(label, tf.int32), depth=3, dtype=logit.dtype)

    loss = tf.nn.softmax_cross_entropy_with_logits(label, logit)

    loss *= mask

    return tf.reduce_sum(loss) / tf.reduce_sum(mask)
