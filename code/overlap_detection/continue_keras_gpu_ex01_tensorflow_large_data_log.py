# -*- coding: utf-8 -*-
"""
Resume training an overlap detection network

Original code by Mighty, repurposed and edited by mkunes
"""

from keras.models import Sequential
from keras.layers.convolutional import Conv2D, MaxPooling2D
from keras.layers.core import Dense, Dropout, Flatten
from keras.layers.normalization import BatchNormalization
from keras.optimizers import SGD, RMSprop
from keras.models import load_model
from keras.callbacks import CSVLogger
import argparse
import h5py
import numpy
import sys
import time
import os
import random
import copy
          
# parse commandline
parser = argparse.ArgumentParser(description='Read data from overlap detection dataset. [\'spektro\'] is a spectrogram (512, 94) dataset and [\'label\'] contains labels with a value between 0 and 1')
parser.add_argument('h5_dir_path', type=str, help='path to a directory with h5 files containing training data, \'spektro\' and \'label\' keys')
parser.add_argument('batch_size', type=int, help='number of samples per batch')
parser.add_argument('max_epoch', type=int, help='number of maximum epoch')
parser.add_argument('update_lr', type=bool, help='True/False - updates lr (*0.1) when loss is stagnant (if SGD)')
parser.add_argument('input_net_name', type=str, help='name of the input network')
parser.add_argument('net_name', type=str, help='name of the output network')
args=parser.parse_args()

# number of images per batch (needs to be round)
batch_size = args.batch_size
# maximal iterations (one iteration = one batch process)
max_iter = args.max_epoch

model = Sequential()

# CNN for spectrogram

model = load_model(args.input_net_name)

sgd = SGD(lr=0.005, decay=0.0005, momentum=0.9, nesterov=True)
#rmsProp = RMSprop(lr=0.00001, rho=0.9, epsilon=1e-06)
#model.compile(loss='binary_crossentropy', optimizer=rmsProp)
model.compile(loss='binary_crossentropy', optimizer=sgd)

print(model.summary())

csv_logger = CSVLogger(args.net_name + '_log.csv', append=True, separator=';')
text_file = open(args.net_name + '_filenames.txt', 'w')

json_string = model.to_json()
open(args.net_name + '_architecture.json', 'w').write(json_string)

h5files = os.listdir(args.h5_file_path)
f_list = []
for h5f in h5files:
    if h5f.endswith('.h5'):
        f_list.append(h5f)

usable_h5 = copy.deepcopy(f_list)

for i in range(args.max_epoch):
    if len(usable_h5) == 0:
       usable_h5 = copy.deepcopy(f_list)

    r = random.randint(0, len(usable_h5)-1)
    h5 = usable_h5[r]

    print(h5)

    f = h5py.File(os.path.join(args.h5_file_path, h5), 'r')
    dX = f['spektro'][()]
    dY = f['label'][()]

    dX = numpy.expand_dims(dX, 1)

    print('Data loaded...')
    print('Training net...')
    print('Epoch: ' + str(i) + '/' + str(args.max_epoch))
    
    model.fit(dX, dY, batch_size=args.batch_size, epochs=1, shuffle=True,validation_split=0.1, callbacks=[csv_logger])
    text_file.write("%s\n" % h5)
    
    del usable_h5[r]
    
    if i % 6 == 0 and i != args.max_epoch and i != 0:
        model.save(args.net_name + '_epoch_' +  str(i) + '.h5')

model.save_weights(args.net_name + '_weights.h5', overwrite=True)
model.save(args.net_name)
text_file.close()
