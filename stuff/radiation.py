import numpy as np
import math
import matplotlib.pyplot as plt

# global constants
c = 2.998*10**8         # Lichtgeschindigkeit in m/s
h = 6.626*10**-34       # Plancks Konstante in Js
k = 1.381*10**-23       # Boltzmankonstante in J/K
sigma = 5.67 *10**-8    # Stefan-Boltzmann in W/ m^2 K^4
S = 1370                # Solarkonstante in W/m^2
AU = 150*10**11         # Astronomische Einheit in m

def micro2m(wav):
    wav_neu = wav*10**-6
    return wav_neu

def freq2wav(freq):
    """frequen in pro sekunde, GHz = 10**9"""
    wave = c / freq
    return wave

def wav2freq(wav):
    """wav in meter"""
    frequenc = c / wav
    return frequenc

def K2C(Kelvin):
    C = Kelvin - 273.15
    return C

def C2K(C):
    Kelvin = C + 273.15
    return Kelvin

def planck(wav, T):
    a1 = (2.0*h*c**2)/(wav**5)
    b1 = (h*c)/(k*T*wav)
    intensity = a1 * 1/(np.exp(b1)-1.0)
    return intensity

def emission(intens, emiss):
    Intens_neu = emiss*intens
    return Intens_neu

def intens2Tb(wav,intens):
    a1 = (2.0*h*c**2)/(intens*wav**5)
    b1 = (h*c)/(k*wav)
    Tb = b1 * 1/(math.log(1.0 + a1))
    return Tb


def wtr(A, T):
    """WattThermalRadiation"""
    W = A * sigma * T**4
    return W

def nhl(A, T1, T2):
    """net Heat loss"""
    W = A * sigma * (T1**4 - T2**4)
    return W

def Teffplanet(A, D):
    Teff = ((S * (1.-A)) / (4.*sigma*D**2.))**(1./4.)
    return Teff



#G.P.?
#plt.plot(w,planck(w,3000))

#GP 6.15
#a) wtr(2,305.2)*0.97
#b) nhl(2, C2K(32), C2K(20))*0.97



# GP 6.11
#T1 = intens2Tb(micro2m(11), planck(micro2m(11),300)*0.95)
#T2 = intens2Tb(0.0158, planck(0.0158,300)*0.95)

# GP 6.20
Ai = np.array([0.11, 0.65, 0.30, 0.15, 0.52, 0.41, 0.30])
Di = np.array([0.39, 0.72, 1.0, 1.52, 5.20, 30.0, 40.0])
for i in range(len(Ai)):
    print Teffplanet(Ai[i] , Di[i])


emi = np.arange(0.,1,0.01)
tempi = np.arange(250,350)
wavi = np.arange(micro2m(0.1),micro2m(7), micro2m(0.01))

for i in emi:
    plt.scatter(i,intens2Tb(0.0158, planck(0.0158,300)*i))


plt.xlabel('Emmision')
plt.ylabel('Temperatur in K')
plt.show()

