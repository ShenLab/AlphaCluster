import tensorflow as tf
import sys
path = sys.argv[1]

model = tf.keras.models.load_model(path)

model.summary()
