
import sys
import wave

import numpy                as np
import matplotlib.pyplot    as plt
import scipy.io.wavfile     as spw


class FrequencyEstimator():
    def __init__(self, frames, frameRate, verbose=False):
        self.verbose = verbose
        self.Fmin = 50
        self.Fmax = 900
        self.H = 4
        self.frameRate = frameRate

        self.N = int(0.7 * frameRate)
        self.Nfft = self.nextPower(self.N, 2)
        self.dFmin = frameRate / self.N

        offset = int(0.1 * frameRate)
        self.x = np.array(frames[offset:offset + self.N]) * np.hamming(self.N)
        self.fft = np.fft.fft(self.x, n=self.Nfft)

    def nextPower(self, n, p):
        '''Outputs the next power of p for n.'''
        m = 1
        while m <= n:
            m *= p
        return m

    def computeProduct(self):
        '''Computes the spectral product.'''

        self.R = int(self.Nfft / 2 / self.H + 1)
        self.p = np.ones(self.R)

        if self.verbose:
            print('R: {}'.format(str(self.R)))

        for r in range(self.R):
            for h in range(self.H):
                self.p[r] *= abs(self.x[r*h])

    def findMax(self):
        '''Computes the maximal frequency of the spectral product.'''

        Nmin = int(self.Fmin / self.frameRate * self.Nfft)
        Nmax = min(self.R, int(self.Fmax / self.frameRate * self.Nfft))

        if self.verbose:
            print('Nmin: {}'.format(str(Nmin)))
            print('Nmax: {}'.format(str(Nmax)))

        i = Nmin + np.argmax(self.p[Nmin:Nmax])
        f = i / self.Nfft

        return f


if __name__ == '__main__':

    with wave.open(sys.argv[1]) as s:
        frames = []
        frame = s.readframes(1)
        while frame:
            i = int.from_bytes(frame, byteorder='big')
            i = int.from_bytes(frame, byteorder='little')
            frames.append(i)
            frame = s.readframes(1)


    print("with scipy")
    rate, frames2 = spw.read(sys.argv[1])
    frameRate = s.getframerate()
    frequencyEstimator = FrequencyEstimator(frames2, frameRate, verbose=True)

    frequencyEstimator.computeProduct()
    f0 = frequencyEstimator.findMax()
    print('Frequence fondamentale : {}'.format(str(f0)))

    plt.figure() 
    plt.subplot(121) ; plt.plot(frames)
    plt.subplot(122) ; plt.plot(frames2)
    plt.tight_layout() ; plt.show()

    print("homecooked")
    frameRate = s.getframerate()
    frequencyEstimator = FrequencyEstimator(frames2, frameRate, verbose=True)

    frequencyEstimator.computeProduct()
    f0 = frequencyEstimator.findMax()
    print('Frequence fondamentale : {}'.format(str(f0)))
