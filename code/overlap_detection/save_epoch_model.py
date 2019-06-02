

from keras.models import Sequential
from keras.layers.convolutional import Conv2D, MaxPooling2D
from keras.layers.core import Dense, Dropout, Flatten
from keras.layers.normalization import BatchNormalization
from keras.optimizers import SGD, RMSprop
from keras.models import load_model
import argparse
import h5py
import numpy
import sys
import time
import os
          
# parse commandline
parser = argparse.ArgumentParser(description='Save an in-progress NN as weights')
parser.add_argument('input_net_name', type=str, help='name of the input network')
parser.add_argument('net_name', type=str, help='name of the output network')
args=parser.parse_args()


model = Sequential()


model = load_model(args.input_net_name)

print(model.summary())

json_string = model.to_json()
open(args.net_name + '_architecture.json', 'w').write(json_string)

model.save_weights(args.net_name + '_weights.h5', overwrite=True)
model.save(args.net_name)


