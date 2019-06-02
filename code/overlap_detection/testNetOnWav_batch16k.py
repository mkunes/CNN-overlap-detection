# -*- coding: utf-8 -*-
"""
Perform overlap detection on 16kHz wav files

Original code by Mighty, repurposed and edited by mkunes
"""
import argparse
from keras.models import model_from_json
from keras.optimizers import RMSprop, SGD

import sys
import os
import numpy
import computeSG
from scipy.io.wavfile import read

# parse commandline
parser = argparse.ArgumentParser(description='Test speaker overlap detection on a wav (16kHz).')
parser.add_argument('wav_list', type=str, help='path to text file with paths to testing waves')
parser.add_argument('input_net_name', type=str, help='name of the input network architecture (json, yaml)')
parser.add_argument('input_net_weights', type=str, help='name of the input network weights (h5)')
parser.add_argument('output', type=str, help='path to output')
args=parser.parse_args()

# load the CNN (as graph)
model = model_from_json(open(args.input_net_name).read())
model.load_weights(args.input_net_weights)

#rmsProp = RMSprop(lr=0.001, rho=0.9, epsilon=1e-06)
sgd = SGD(lr=0.01, decay=0.0005, momentum=0.9, nesterov=True)
model.compile(loss='binary_crossentropy', optimizer=sgd)

print(model.summary())

f = open(args.wav_list, 'r')
waves = f.readlines()
f.close()

delta_t = 1 #1 #1.4
shift_t = 0.05
freq = 1024 # use 512 for 8kHz, 1024 for 16kHz

for w in waves:
    wave_info = read(w.rstrip())
    wave = wave_info[1].astype(numpy.float32) / 32768

    if len(wave_info[1].shape) > 1:
        if wave_info[1].shape[1] > 1:
            wave = wave_info[1][:, 0].astype(numpy.float32) / 32768

    Fs = wave_info[0]

    result = []

    progress = 0.0

    t = 0.0

    p, f = os.path.split(w)
    fn, fe = os.path.splitext(f)
    fn = fn + '.bin'

    while True:
        sg = numpy.log(computeSG.computeSG(wave[int(numpy.floor(t*Fs)):int(numpy.floor((t+delta_t)*Fs))], freq, int(numpy.floor(Fs/100))))
        #sg = sg[:, 0:256]
        sg = numpy.transpose(sg)

        val = model.predict(numpy.array([[sg]]), batch_size=1, verbose=0)
        #print str(val)
        result.append(val)

        t = t + shift_t
        if (t + delta_t)*Fs > wave.shape[0]:
            break

        # progress = (t + delta_t)*Fs / wave.shape[0] * 100.0
        # print ('progress: ' + str(progress) + ' precent\r')
        #
        # sys.stdout.flush()

    numpy.array(result, dtype=numpy.float32).tofile(os.path.join(args.output, fn))