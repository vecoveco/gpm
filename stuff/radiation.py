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
R_Erde = 6356.0         # Erdradius in km
kw = 2897.0             # Wien's Komstante in nm K


def BeerBougetLambert(F, beta, s):
    """
    Function:
         bouguer-lambert-beersche Gesetz beschreibt die Abschwaechung der
         Intensitaet einer Strahlung bei dem Durchgang durch ein Medium
         mit einer absorbierenden Substanz, in Abhaengigkeit von der
         Konzentration der absorbierenden Substanz und der Schichtdicke
    Input:
        F       ::: radiation flux density in W/m2 befor absorption
        beta    ::: volumen absorption coeff.
        s       ::: path length
    Output:
        F_att   ::: radiation flux density in W/m2 after absorption

    """
    F_att = F * np.exp(-1* beta * s)
    return F_att


def beta(ni, lam):
    """
    Function:
         Berechnung des volumen absorptions coeff
    Input:

        ni    ::: refractive Index
        lam   ::: wellenlaenge
    Output:
        beta  ::: absorptions coeff

    """
    beta = (4* math.pi * ni)/lam
    return beta

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

def Teff(A, s):
    Teff = ((s * (1.-A)) / (4.*sigma))**(1./4.)
    return Teff

def SB(S):
    T = (S/sigma)**(1./4.)
    return T

def SB2(T):
    S = sigma*T**4.
    return S

def srtm(asw, alw, A, eps=None):
    """
    Simple radiative transfere model of Atmo
    A ::: Albedo
    asw ::: Absorption kurzwellig
    alw ::: Absorbtion langwellig
    """
    if eps == None:
        Sm = S / 4.
        Esurf = ((1. - (1. - asw) * A) * ((2. - asw)/(2. - alw)))
        Tsurf = (Sm / sigma * (Esurf))**(1./4.)
        Eatmo = (((1 - A) * (1 - asw) * alw) + (1 + (1 - asw) * A) * asw)/((2 - alw) * alw)
        Tatmo = (Sm/sigma * (Eatmo))**(1./4.)
        return Tatmo, Tsurf
    if eps != None:
        eps = eps
        Sm = S / 4.
        Esurf = ((A - 1) * (1 - asw) * (1 + (1 - alw)) * (1 - eps) + eps * (A * (1 - asw)**2)-1)/(eps * (2 - alw))
        Tsurf = (-1* Sm / sigma * (Esurf))**(1./4.)
        Eatmo = (((A - 1) * (1 - asw))*alw + (1 + ((1 - asw) * A))* asw) / (alw * ( eps * (( 1- alw) + 1)) + ((1 - alw) * (1-eps)))
        Tatmo = (-1*Sm/sigma * (Eatmo))**(1./4.)
        return Tatmo, Tsurf



#Ex6.19a
a_sw = np.arange(0.01,1.0,0.01)
a_lw = np.arange(0.01,1.0,0.01)
Albedo = np.arange(0.01, 1.0, 0.01)
emi = np.arange(0.01, 1.0, 0.01)


# Ex 19a
Ta, Ts = srtm(a_sw, 0.8, 0.3)
plt.plot(a_sw, Ta, label='Ta', color='blue')
plt.plot(a_sw, Ts, label='Ts', color='green' )
plt.ylabel('Temp in K')
plt.xlabel('Absorbiton')
plt.grid()
plt.legend(loc='upper right')
plt.show()

# Ex 19b
Ta, Ts = srtm(0.1, a_lw, 0.3)
plt.plot(a_sw, Ta, label='Ta', color='blue')
plt.plot(a_sw, Ts, label='Ts', color='green' )
plt.ylabel('Temp in K')
plt.xlabel('Absorbiton')
plt.grid()
plt.legend(loc='upper right')
plt.show()

# Ex 19c
Ta, Ts = srtm(0.1, 0.8, Albedo)
plt.plot(a_sw, Ta, label='Ta', color='blue')
plt.plot(a_sw, Ts, label='Ts', color='green' )
plt.ylabel('Temp in K')
plt.xlabel('Albedo')
plt.grid()
plt.legend(loc='upper right')
plt.show()

for i in Albedo:
    Ta, Ts = srtm(a_sw, 0.8, i)
    plt.plot(a_sw, Ta, label='Ta', color='blue')
    plt.plot(a_sw, Ts, label='Ts', color='green' )


plt.ylabel('Temp in K')
plt.xlabel('Absorbiton')
plt.grid()
plt.legend(loc='upper right')
plt.show()

# Ex 19d
Ta, Ts = srtm(0.1, 0.8, 0.3, emi)
plt.plot(emi, Ta, label='Ta', color='blue')
plt.plot(emi, Ts, label='Ts', color='green' )
plt.ylabel('Temp in K')
plt.xlabel('Albedo')
plt.grid()
plt.legend(loc='upper right')
plt.show()


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

