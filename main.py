from scipy.io.wavfile import read
from scipy.io.wavfile import write
from math import log10
from math import cos
from math import pi
import matplotlib.pyplot as plt
import numpy as np
import numpy.fft as FFT

def main(arg):
    rate, data = read(arg)
    affiche(rate, data)

def affiche(rate, data):
    print("frequence d'echantillonnage : {0} Hz\n".format(rate))
    print("taille fichier :\n{0} echantillons".format(len(data)))
    print("{0} ms\n".format(len(data)*1000/rate))
    plt.subplot(311)
    plt.plot(data)
    signal_modif = modif_signal(rate, data, 8*rate/1000, 32*rate/1000)
    write("resultat.wav", rate, np.int16(signal_modif))
    plt.subplot(312)
    plt.plot(signal_modif)

    plt.show()

def fenetrageHamming(size):
    hamming = []
    for n in range(0, size):
        angle = 2 * pi * n/size
        hamming.append(0.54 - 0.46 * cos(angle))

    return hamming

def fenetrage(signal, hamming):
    fen = []
    for i in range(0, len(hamming)):
        fen.append(signal[i]*hamming[i])

    return fen

def spectre_amplitude(spectre, fftsize):
    res = np.abs(spectre)
    spec_affichage = res
    res = [20*x for x in res]
    spec_affichage = [20*log10(x) for x in spec_affichage]

    return res, spec_affichage
    
def modif_signal(rate, signal, m, N):
    signal_modif = np.zeros(len(signal))
    somme_hamming = np.zeros(len(signal))
    tab_spectres = [] # pour le stockage des spectres d'amplitude
    hamming = fenetrageHamming(N)
    fftsize = 1024
    
    for i in range(0, len(signal) - N, m):
        fenetre = np.array(fenetrage(signal[i:i+N], hamming), dtype=np.float)
        spectre = FFT.fft(fenetre, fftsize)
        amplitude, ampli_aff = spectre_amplitude(spectre, fftsize)
        tab_spectres.append(ampli_aff)
        fenetre = np.real(FFT.ifft(spectre, fftsize))
        signal_modif[i:i+N] += fenetre[0:N]
        somme_hamming[i:i+N] += hamming

    np.transpose(tab_spectres)
    plt.subplot(313)
    plt.imshow(tab_spectres, aspect='auto') # afficher uniquement 0 a fftsize/2
    signal_modif /= somme_hamming
    return signal_modif
        
    
if __name__ == "__main__":
    main("test_seg.wav")
