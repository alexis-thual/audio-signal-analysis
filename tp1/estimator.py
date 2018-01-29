
import numpy as np
import matplotlib.pyplot as plt


class FrequencyEstimator:

    def __init__(self, frames, frameRate, verbose=False):
        self.verbose = verbose
        self.colors = ["#2965CC", "#29A634", "#D99E0B", "#D13913"]

        self.Fmin = 50
        self.Fmax = 900
        self.H = 4
        self.frameRate = frameRate

        self.N = int(0.7 * frameRate)
        self.Nfft = self.nextPower(self.N, 2)
        self.dFmin = frameRate / self.N

        offset = int(0.1 * frameRate)
        self.x = np.array(frames[offset:offset + self.N]) * np.hamming(self.N)
        self.x /= abs(self.x.max() + 1e-10)
        self.fft = np.fft.fft(self.x, n=self.Nfft)
        self.fft = np.absolute(self.fft)

        self.alpha = 0.02
        self.beta = 0

        if self.verbose:
            self.plotSignal()

    def plotSignal(self):
        f, axarr = plt.subplots(1, 2, figsize=(15, 5))
        axarr[0].plot(self.x, color=self.colors[0])
        axarr[0].set_title("X through Hamming fct")
        axarr[1].plot(self.fft[10:], color=self.colors[0])
        axarr[1].plot(np.mean(self.fft), color=self.colors[1])
        axarr[1].set_yscale('log')
        axarr[1].set_title("Through Fourier")
        plt.show()

    def nextPower(self, n, p):
        '''Outputs the next power of p for n.'''
        m = 1
        while m <= n: m *= p
        return m

    def computeProduct(self):
        '''Computes the spectral product.'''

        self.R = int(self.Nfft / 2 / self.H + 1)
        self.p = np.ones(self.R)

        if self.verbose:
            print('R: {}'.format(str(self.R)))

        for r in range(self.R):
            for h in range(self.H):
                self.p[r] *= self.fft[min(r * 2**h, self.R - 1)]

        if self.verbose:
            plt.figure(figsize=(15, 5))
            plt.plot(self.p[1:], color=self.colors[0])
            plt.title("Product")
            plt.yscale('log')
            plt.show()

    def findMax(self):
        '''Computes the maximal frequency of the spectral product.'''

        Nmin = int(self.Fmin / self.frameRate * self.Nfft)
        Nmax = min(self.R, int(self.Fmax / self.frameRate * self.Nfft))

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

        k1 = int(f1 * self.Nfft / self.frameRate)
        k2 = int(f2 * self.Nfft / self.frameRate)

        sloc = k1 + np.argmax(self.p[k1:k2])
        swin = 25

        if self.verbose:
            print('f1 {} f2 {}\nk1 {} k2 {}\nspikeLocation {}'.format(str(f1), str(f2), str(k1), str(k2), str(sloc)))
            print('Corrected frequency: {}'.format(str(sloc * self.frameRate / self.N)))

        kf = sloc
        while kf < self.Nfft:
            self.fft[sloc - swin:sloc + swin] = np.repeat(self.fft[sloc - swin:sloc + swin].min(), 2 * swin)
            kf += sloc

    def shouldStop(self):
        return False
