from scipy.io.wavfile import read
from scipy.io.wavfile import write
import math
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

def fenetrage_hamming(size):
    hamming = []
    for n in range(0, size):
        angle = 2 * math.pi * n/size
        hamming.append(0.54 - 0.46 * math.cos(angle))

    return hamming

def fenetrage(signal, hamming):
    fen = []
    for i in range(0, len(hamming)):
        fen.append(signal[i]*hamming[i])

    return fen

def spectre_amplitude(spectre):
    res = np.abs(spectre)
    spec_affichage = res
    spec_affichage = [20*math.log10(x) for x in spec_affichage]

    return res, spec_affichage

def spectre_phase(spectre):
    return np.angle(spectre)

def spectre_reconstruction(spectre_amplitude, spectre_phase):
    return spectre_amplitude * np.exp(1j * spectre_phase)

def soustraction_spectrale(spectre_amplitude, bruit):
    alpha = 2
    beta = 8
    gamma = 0.2
    res = np.zeros(len(spectre_amplitude))
    y = beta * math.pow(bruit, alpha)
    for i in range(0, len(spectre_amplitude)):
        x = math.pow(spectre_amplitude[i], alpha)
        if (x - y) > 0:
            res[i] = math.pow(x - y, 1.0/alpha)
        else:
            res[i] = gamma * bruit

    return res
        
def modif_signal(rate, signal, m, N):
    signal_modif = np.zeros(len(signal))
    somme_hamming = np.zeros(len(signal))
    tab_spectres = [] # pour le stockage des spectres d'amplitude
    hamming = fenetrage_hamming(N)
    fftsize = 1024

    # estimation du bruit
    bruit = 0
    size = 5
    for i in range(0, (size - 1)*m, m):
        fenetre = np.array(fenetrage(signal[i:i+N], hamming), dtype=np.float)
        spectre = FFT.fft(fenetre, fftsize)
        amplitude, ampli_aff = spectre_amplitude(spectre)
        bruit += np.mean(amplitude)
    bruit /= size
    
    for i in range(0, len(signal) - N, m):
        fenetre = np.array(fenetrage(signal[i:i+N], hamming), dtype=np.float)
        spectre = FFT.fft(fenetre, fftsize)
        amplitude, ampli_aff = spectre_amplitude(spectre)
        tab_spectres.append(ampli_aff)
        phase = spectre_phase(spectre)
        # modification spectre d'amplitude
        amplitude = soustraction_spectrale(amplitude, bruit)
        # reconstruction du signal
        spectre = spectre_reconstruction(amplitude, phase)
        fenetre = np.real(FFT.ifft(spectre, fftsize))
        signal_modif[i:i+N] += fenetre[0:N]
        somme_hamming[i:i+N] += hamming

    np.transpose(tab_spectres)
    plt.subplot(313)
    plt.imshow(tab_spectres, aspect='auto') # affichage spectrogramme
    signal_modif /= somme_hamming
    return signal_modif
        
    
if __name__ == "__main__":
    main("test_seg_bruit_0dB.wav")
