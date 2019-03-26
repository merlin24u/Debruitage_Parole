from scipy.io.wavfile import read
from scipy.io.wavfile import write
from math import log10
from math import cos
from math import pi
import matplotlib.pyplot as plt
import numpy as np

def main(arg):
    rate, data = read(arg)
    affiche(rate, data)

def affiche(rate, data):
    print("frequence d'echantillonnage : {0} Hz\n".format(rate))
    print("taille fichier :\n{0} echantillons".format(len(data)))
    print("{0} ms\n".format(len(data)*1000/rate))
    plt.subplot(311)
    plt.plot(data)
    hamming = fenetrageHamming(32)
    fenetrage(data, hamming)
    plt.show()

def fenetrageHamming(size):
    hamming = []
    for n in range(0, size):
        angle = 2 * pi * n/size
        hamming.append(0.54 - 0.46 * cos(angle))
    
    plt.subplot(312)
    plt.plot(hamming)

    return hamming

def fenetrage(signal, hamming):
    fen = []
    for i in range(0,len(hamming)):
        fen.append(signal[i]*hamming[i])

    plt.subplot(313)
    plt.plot(fen)

    return fen

def modif_signal()
    
if __name__ == "__main__":
    main("test_seg.wav")
