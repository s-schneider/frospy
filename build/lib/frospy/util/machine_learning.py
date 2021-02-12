import tensorflow as tf
from tensorflow.keras import layers
from frospy.util.read import read_pickle
from os.path import join
import numpy as np


def read_data_labels(file, Nsample, label, spectrum='abs'):
    data = read_pickle(file)
    for i, d in enumerate(data):
        if spectrum == 'abs':
            data[i] = abs(d[0:Nsample])
        elif spectrum == 'real':
            data[i] = [x.real for x in d[0:Nsample]]
        elif spectrum == 'imag':
            data[i] = [x.imag for x in d[0:Nsample]]
        else:
            raise IOError('no spectrum defined')
    data = np.array(data).reshape(len(data), Nsample)
    if label == 1:
        labels = np.ones(len(data), dtype='int').reshape(len(data), 1)
    else:
        labels = np.zeros(len(data), dtype='int').reshape(len(data), 1)
    return data, labels


model = tf.keras.Sequential()
# Adds a densely-connected layer with 64 units to the model:
model.add(layers.Dense(64, activation='relu'))
# Add another:
model.add(layers.Dense(64, activation='relu'))
# Add a softmax layer with 10 output units:
model.add(layers.Dense(2, activation='softmax'))

# layers.Dense(64, activation='sigmoid')
model.compile(optimizer=tf.train.AdamOptimizer(0.001),
              loss='categorical_crossentropy',
              metrics=['accuracy'])

# Data Feed
# import numpy as np
#
# data = np.random.random((1000, 32))
# labels = np.random.random((1000, 10))
data_dir = '/quanta1/home/simons/splitting/alldatafiles/tensorflow'
fh_val = 'valid/VHZ.030385A.spectra.pickle'
fh_inval = 'invalid/VHT.010196C.spectra.pickle'
Nsample = 14

file = join(data_dir, fh_val)
data1, labels1 = read_data_labels(file, Nsample, 1)

file = join(data_dir, fh_inval)
data2, labels2 = read_data_labels(file, Nsample, 0)

data = np.vstack((data1, data2))
labels = np.vstack((labels1, labels2))


model.fit(data, labels, epochs=10, batch_size=32)
