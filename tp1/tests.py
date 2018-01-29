
import sys

import matplotlib.pyplot    as plt
import scipy.io.wavfile     as spw

import estimator


rate, frames = spw.read(sys.argv[1])
fe = estimator.FrequencyEstimator(frames, rate, verbose=True)

fe.computeProduct()
f0 = fe.findMax()
print('Frequence fondamentale : {}'.format(str(f0)))
print("Elimination")
fe.eliminateFrequency(f0)
fe.computeProduct()
