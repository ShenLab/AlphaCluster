import tensorflow as tf


class ModelBase(tf.keras.Model):
    def __init__(self):
        super(ModelBase, self).__init__()

    def call(self, inputs, training):
        raise NotImplementedError('The call method has to be override')

    def predict_from_logit(self, logit):
        return tf.nn.softmax(logit)
