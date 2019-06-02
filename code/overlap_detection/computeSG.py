# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 11:33:43 2016

@author: Mighty
"""

import numpy

def computeSG(wav, lWindow, shift):
    half_win = int(numpy.floor(lWindow/2))
    
    hammingWindow = 0.54 - 0.46*numpy.cos(2*numpy.pi*numpy.linspace(0, (lWindow - 1) / (lWindow - 1), lWindow))
    
    Lw = len(wav)
    
    i = 0
    
    sg = numpy.zeros((int(((wav.shape[0]-lWindow)//shift)+1), half_win), dtype=numpy.float32)
    sg_idx = 0
    
    while i + lWindow < Lw:
        win = wav[i:i+lWindow]
        s = numpy.fft.fft(win * hammingWindow)
        
        s2 = s[0:half_win]
        s2 = numpy.abs(s2)

        sg[sg_idx] = s2
        sg_idx += 1
        
        i = i + shift
        
    sg = numpy.reshape(sg, (-1, half_win))
    return sg

def computeLabels(ann, samples, lWindow=512, shift=80, Fs=8000, tau=0.6):
    offset_in_s = lWindow / Fs / 2
    shift_in_s = shift/Fs
    Lw = len(samples)

    simple_ann = numpy.zeros((len(ann)*2))
    i = 0
    for a in ann:
        simple_ann[i] = a[0]
        i += 1
        simple_ann[i] = a[1]
        i += 1

    i = 0
    labels_idx = 0
    labels = numpy.zeros((int(((Lw-lWindow)//shift)+1)), dtype=numpy.float32)

    while i + lWindow < Lw:
        t = labels_idx*shift_in_s + offset_in_s
        v = numpy.min(numpy.abs(t - simple_ann))
        labels[labels_idx] = max(-0.5*v/(tau/2) + 1, 0)

        labels_idx += 1
        i += shift

    return labels

def computeLabelsMFCC(ann, samples, lWindow=0.032, shift=0.016, tau=0.6):
    offset_in_s = lWindow / 2
    Lw = samples.shape[0]

    simple_ann = numpy.zeros((len(ann)*2))
    i = 0
    for a in ann:
        simple_ann[i] = a[0]
        i += 1
        simple_ann[i] = a[1]
        i += 1

    i = 0
    labels_idx = 0
    labels = numpy.zeros((samples.shape[0]), dtype=numpy.float32)

    for i in range(Lw):
        t = labels_idx*shift + offset_in_s
        v = numpy.min(numpy.abs(t - simple_ann))
        labels[labels_idx] = max(-0.5*v/(tau/2) + 1, 0)

        labels_idx += 1

    return labels
    
def computeSG_withPhase(wav, lWindow, shift):
    half_win = int(numpy.floor(lWindow/2))
    
    hammingWindow = 0.54 - 0.46*numpy.cos(2*numpy.pi*numpy.linspace(0, (lWindow - 1) / (lWindow - 1), lWindow))
    
    Lw = len(wav)
    
    i = 0
    
    sgExists = False
    
    while i + lWindow < Lw:
        win = wav[i:i+lWindow]
        s = numpy.fft.fft(win * hammingWindow)
        
        s2 = s[0:half_win]        

        if sgExists == False:
            sg = numpy.zeros((2, half_win, 1))
            sg[0, :, 0] = s2.real
            sg[1, :, 0] = s2.imag
            sgExists = True
        else: 
            s_real = s2.real.reshape(1, -1, 1)
            s_imag = s2.imag.reshape(1, -1, 1)
            
            sc = numpy.concatenate([s_real, s_imag], axis=0)                        
            sg = numpy.concatenate([sg, sc], axis=2)
        
        i = i + shift
        
    #sg = numpy.reshape(sg, (-1, lWindow/2))
    return sg
    
def computeSG_withPhase2(wav, lWindow, shift):
    half_win = int(numpy.floor(lWindow/2))
    
    hammingWindow = 0.54 - 0.46*numpy.cos(2*numpy.pi*numpy.linspace(0, (lWindow - 1) / (lWindow - 1), lWindow))
    
    Lw = len(wav)
    
    i = 0
    
    sgExists = False
    
    while i + lWindow < Lw:
        win = wav[i:i+lWindow]
        s = numpy.fft.fft(win * hammingWindow)
        
        s2 = s[0:half_win]        
        sMag = numpy.abs(s2)

        if sgExists == False:
            sg = numpy.zeros((3, lWindow/2, 1))            
                        
            sg[0, :, 0] = sMag
            sg[1, :, 0] = s2.real / sMag
            sg[2, :, 0] = s2.imag / sMag
            sgExists = True
        else:             
            s_real = s2.real /  sMag;
            s_real = s_real.reshape(1, -1, 1)
            
            s_imag = s2.imag / sMag
            s_imag = s_imag.reshape(1, -1, 1)
            
            s_Mag = sMag.reshape(1, -1, 1)                        
            
            sc = numpy.concatenate([s_Mag, s_real, s_imag], axis=0)                        
            sg = numpy.concatenate([sg, sc], axis=2)
        
        i = i + shift
        
    #sg = numpy.reshape(sg, (-1, lWindow/2))
    return sg
    
def computeSG_withPhase3(wav, lWindow, shift):
    half_win = int(numpy.floor(lWindow/2))
    
    hammingWindow = 0.54 - 0.46*numpy.cos(2*numpy.pi*numpy.linspace(0, (lWindow - 1) / (lWindow - 1), lWindow))
    
    Lw = len(wav)
    
    i = 0
    
    sgExists = False
    
    while i + lWindow < Lw:
        win = wav[i:i+lWindow]
        s = numpy.fft.fft(win * hammingWindow)
        
        s2 = s[0:(lWindow/2)]        
        sMag = numpy.abs(s2)

        if sgExists == False:
            sg = numpy.zeros((3, half_win, 1))            
                        
            sg[0, :, 0] = numpy.log(sMag)
            sg[1, :, 0] = s2.real / sMag
            sg[2, :, 0] = s2.imag / sMag
            sgExists = True
        else:             
            s_real = s2.real /  sMag;
            s_real = s_real.reshape(1, -1, 1)
            
            s_imag = s2.imag / sMag
            s_imag = s_imag.reshape(1, -1, 1)
            
            s_Mag = numpy.log(sMag.reshape(1, -1, 1))
            
            sc = numpy.concatenate([s_Mag, s_real, s_imag], axis=0)                        
            sg = numpy.concatenate([sg, sc], axis=2)
        
        i = i + shift
        
    #sg = numpy.reshape(sg, (-1, lWindow/2))
    return sg