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
    modif_signal(rate, data, 8*rate/1000, 32*rate/1000)
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

def modif_signal(rate, signal, m, N):
    signal_modif = []
    for i in range(0, len(signal) - N, m):
        hamming = fenetrageHamming(N/rate*1000)
        fenetre = fenetrage(signal[i:i+N], hamming)
        
    
if __name__ == "__main__":
    main("test_seg.wav")
