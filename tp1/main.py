import wave
import numpy as np
from math import floor
import matplotlib.pyplot as plt

class FrequencyEstimator():
    def __init__(self, frames, frameRate, verbose=False):
        self.verbose = verbose
        self.colors = ["#2965CC", "#29A634", "#D99E0B", "#D13913"]

        self.Fmin = 50
        self.Fmax = 900
        self.H = 4
        self.frameRate = frameRate

        self.N = int(floor(0.7 * frameRate))
        self.Nfft = self.nextPower(self.N, 2)
        self.dFmin = frameRate / self.N

        offset = int(floor(0.1 * frameRate))
        self.x = np.array(frames[offset:offset + self.N]) * np.hamming(self.N)
        self.fft = np.fft.fft(self.x, n=self.Nfft)
        self.fft = np.absolute(self.fft)

        self.alpha = 0.02
        self.beta = 0

        if self.verbose:
            self.plotSignal()

    def plotSignal(self):
        f, axarr = plt.subplots(1, 2, figsize=(15,5))
        axarr[0].plot(self.x, color=self.colors[0])
        axarr[1].plot(self.fft[10:], color=self.colors[0])
        axarr[1].plot(np.mean(self.fft), color=self.colors[1])
        axarr[1].set_yscale('log')
        plt.show()

    def nextPower(self, n, p):
        '''Outputs the next power of p for n.'''
        m = 1
        while m <= n:
            m *= p
        return m

    def computeProduct(self):
        '''Computes the spectral product.'''

        self.R = int(floor(self.Nfft / 2 / self.H + 1))
        self.p = np.ones(self.R)

        if self.verbose:
            print('R: {}'.format(str(self.R)))

        for r in range(self.R):
            for h in range(self.H):
                self.p[r] *= self.fft[r*h]

        if self.verbose:
            fig = plt.figure(figsize=(15,5))
            # Eliminate frequency 0
            plt.plot(self.p[1:], color=self.colors[0])
            plt.yscale('log')
            plt.show()

    def findMax(self):
        '''Computes the maximal frequency of the spectral product.'''

        Nmin = int(floor(self.Fmin / self.frameRate * self.Nfft))
        Nmax = min(self.R, int(floor(self.Fmax / self.frameRate * self.Nfft)))

        if self.verbose:
            print('Nmin: {}'.format(str(Nmin)))
            print('Nmax: {}'.format(str(Nmax)))

        i = Nmin + np.argmax(self.p[Nmin:Nmax])
        f = i * self.frameRate / self.Nfft

        if self.verbose:
            print('Found frequency: {}'.format(str(f)))

        return f

    def eliminateFrequency(self, f):
        f1 = (1 - self.alpha) * f
        f2 = (1 + self.alpha) * f

        k1 = int(floor(f1 * self.Nfft / self.frameRate))
        k2 = int(floor(f2 * self.Nfft / self.frameRate))

        spikeLocation = k1 + np.argmax(self.p[k1:k2])
        spikeSize = 50

        if self.verbose:
            print('f1 {} f2 {} k1 {} k2 {} spikeLocation {}'.format(str(f1), str(f2), str(k1), str(k2), str(spikeLocation)))
            print('Corrected frequency: {}'.format(str(spikeLocation * self.frameRate / self.N)))

        kf = spikeLocation
        while kf < self.Nfft:
            self.fft[spikeLocation - spikeSize/2:spikeLocation + spikeSize/2] = [np.min(self.fft[spikeLocation - spikeSize/2:spikeLocation + spikeSize/2])] * spikeSize
            kf += spikeLocation

    def shouldStop(self):
        return False


if __name__ == '__main__':
    with wave.open('./tp1/sons_multipitch/A3_piano.wav') as s:
        frames = []
        frame = s.readframes(1)
        while frame:
            i = int.from_bytes(frame, byteorder='big')
            i = int.from_bytes(frame, byteorder='little')
            frames.append(i)
            frame = s.readframes(1)

        frameRate = s.getframerate()
        frequencyEstimator = FrequencyEstimator(frames, frameRate, verbose=True)

        frequencies = []
        shouldStop = frequencyEstimator.shouldStop()

        k = 0
        while not shouldStop and k < 10:
            frequencyEstimator.computeProduct()
            f0 = frequencyEstimator.findMax()
            frequencyEstimator.eliminateFrequency(f0)

            frequencies.append(f0)
            shouldStop = frequencyEstimator.shouldStop()
            k += 1

            print('Frequence fondamentale : {}'.format(str(f0)))

        s.close()

# %% Test cell
import numpy as np

a = [1,2,3,4]
a[1:3] = [7,8]
a
