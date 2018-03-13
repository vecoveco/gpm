#!/usr/bin/env python
# -*- coding: UTF-8 -*-


"""

This is a self made Library for the Lecture :

- Introduction to Atmospheric Radiation and Remote Sensing

- The whole contents is from the Lecture 416 and the Book:

> A First Course In Atmospheric Radiation
> Secon edition
> Grant W. Petty


Creat by V. Pejcic

"""

# Import Librarys
import numpy as np
import math
import matplotlib.pyplot as plt


# Important Global Constants
c = 2.998*10**8         # Lichtgeschindigkeit in m/s
h = 6.626*10**-34       # Plancks Konstante in Js
k = 1.381*10**-23       # Boltzmankonstante in J/K
sigma = 5.67 *10**-8    # Stefan-Boltzmann in W/ m^2 K^4
solar = 1370            # Solarkonstante in W/m^2
au = 150*10**11         # Astronomische Einheit in m
r_earth = 6356.0        # Erdradius in km
kw = 2897.0             # Wien's Komstante in nm K
n_a = 6.022*1e23        # Avogadro Number mole^-1


# Frequenz von bestimmten Messgeräten
freq_ku_band = 13.6  #GHz, 1e9 in Hz DPR GPM Ku-band
freq_ka_band = 35.5  #GHz, 1e9 in Hz DPR GPM Ka-band
freq_x_band = 9.3    #GHz, 1e9 in Hz BoXPol/JuXPol Uni Bonn
freq_c_band = 5.6    #GHz, 1e9 in Hz DWD Regenradar
freq_k_band = 24.1   #GHz, 1e9 in Hz # MRR in Bonn
freq_s_band = 3.     #GHz, 1e9 in Hz # Kein Bestimmtes S-band


# Wichtige Funktionen der Atmosphärischen Strahlung
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
    """
    Function:
        Umrechnung von mikrometer in meter.
    Input:
        wav Wellenlaenge in mikrometer
    Output:
        Wellenleange in meter
    """
    wav_neu = wav*10**-6
    return wav_neu

def m2micro(wav):
    """
    Function:
        Umrechnung von meter in mikrometer.
    Input:
        wav Wellenlaenge in meter
    Output:
        Wellenleange in mikrometer
    """
    wav_neu = wav*10**6
    return wav_neu

def freq2wav(freq):
    """
    Funktion: 
        Frequenz in Wellenlaenge umwandeln
        
    Input:
        Frequenz in pro sekunde bzw Herz [Hz] angeben
        
    Output:
        Welenlaenge in Metern [m]
    """
    wave = c / freq
    return wave

def wav2freq(wav):
    """
    Funktion: 
        Wellenlaenge in Frequenz umwandeln
        
    Input:
        Wellenlaenge in Metern
    Output:
        Frequenz in pro sekunde in Herz [Hz]
    """
    frequenc = c / wav
    return frequenc

def K2C(Kelvin):
    """
    Funktion:
        Umrechnung von Kelvin zu Celsius
    Input:
        Temperatur in Kelvin [K]
    Output:
        Temperatur in Celsius [C]
    """
    C = Kelvin - 273.15
    return C

def C2K(C):
    """
    Funktion:
        Umrechnung von Celsius in Kelvin
    Input:
        Temperatur in Celsius [C]
    Output:
        Temperatur in Kevin [K]
    """
    Kelvin = C + 273.15
    return Kelvin

def planck(wav, T):
    """
    Diese Funktion berechnet die Strahlung nach Panlk mit Abhänngigkeit von
    Temperatur und Wellenlänge eines Schwarzen Körpers GP (6.1)
        
    Input:
        wav : Wellenlänge in [µm]
        T : Temperatur in [K]
        
    Output:
        B = Strahlung in [W/(m^2 sr µm)]
    """
    a1 = (2.0*h*c**2)/(wav**5)
    b1 = (h*c)/(k*T*wav)
    B = a1 * 1/(np.exp(b1)-1.0)
    return B

def emission(intens, emiss):
    """
    Monochromatische Emission eines grauen Körpers
        
    Input:
        intens : Strahldichte in [W/m^2]
        emiss : Emissionsgrad in [/]
    Output:
        Intens_neu : Ausgestrahlte Strahldichte in [W/m^2 sr]
        
    """
    Intens_neu = emiss*intens
    return Intens_neu

def intens2Tb(wav,intens):   
    """
    Berechnung der Strahlungstemperatur aus der Plank Funktion.
    GP (6.13)
        
    Input:
        wav : Wellenlänge in [µm]
        intens : Strahldichte in [W/m^2 sr]
    Output:
        Tb : Strahlungstemperatur in [K]
    """
    a1 = (2.0 * h * c**2) / (intens * wav**5)
    b1 = (h * c) / (k * wav)
    Tb = b1 * 1 / (math.log(1.0 + a1))
    return Tb

def stefan_boltzmann(T):
    """
    Stefan-Boltzman_Gesetz Berechnet den Strahlungsflussdichte integriert über
    alle Wellenlängen in abhänngigkeit von der Temperatur T. GP (6.5)
        
    Input:
        T : Temperatur in [K]

    Output:
        Strahlungsflussdichte in [W/m^2]
        
    """
    Fbb = sigma * T**4.
    return Fbb

def srtm(asw, alw, A, eps=None):
    """
    Simpel Radiative Models of Atmosphere zur Berechnung von Atmosphären und
    Oberflächen Temperaturen unter verschiedenen Emissionsgraden. GP S. 140-143
        
    Input:
        A : Albedo in [/]
        asw : Absorption kurzwellig in [/]
        alw : Absorbtion langwellig in [/]
        
    Output:
        Tatmo : Atmosphärentemperatur in [/]
        Tsurf : Oberflächentemperatur in [/]

    """
    if eps == None:
        Sm = S / 4.
        Esurf = ((1. - (1. - asw) * A) * ((2. - asw) / (2. - alw)))
        Tsurf = (Sm / sigma * Esurf)**(1./4.)
        Eatmo = (((1. - A) * (1. - asw) * alw) + (1. + (1. - asw) * A) * asw)\
                / ((2. - alw) * alw)
        Tatmo = (Sm / sigma * Eatmo)**(1./4.)
        return Tatmo, Tsurf
    if eps != None:
        eps = eps
        Sm = S / 4.
        Esurf = ((A - 1.) * (1. - asw) * (1. + (1. - alw)) * (1. - eps)
                 + eps * (A * (1. - asw)**2)-1.) / (eps * (2. - alw))
        Tsurf = (-1. * Sm / sigma * Esurf)**(1. / 4.)

        Eatmo = (((A - 1.) * (1. - asw)) * alw + (1. + ((1. - asw) * A)) * asw)\
                / (alw * (eps * ((1. - alw) + 1.)) + ((1. - alw) * (1. - eps)))

        Tatmo = (-1. *Sm / sigma * Eatmo)**(1. / 4.)
        return Tatmo, Tsurf
    

def rainbow_min_def_ang(m, k):
    """
    Computes a minimum angle of deflection for a Rainbow
    Returns

    Input:
        m : Brechungsindex in [/]
        k : Anzahl der Brechungen
    Output:
        res : Minimaler Winkel der Deflection
    """
    res = np.sqrt((m**2)/(k**2 +2*k))
    res = np.rad2deg(np.arccos(res))

    return res


def optt_cloud(N,H=100, L=0.01, roh=1000, Q=2):
    """
    Funktion zur Berechnung der optischen Dicke einer Wolke
    
    Input: 
        H : Hohe der Wolke in [m]
        N : Anzahl der Tropfen in [1/m^3]
        L : Flüssigwassergehalt in [kg/m^2]
        rho : Dichte in [kg/m^3]
    
    Output: 
        tau : optische Dicke in [/]
    
    """
    ne = 9. * H * np.pi * N * (L**2)  
    za = 16. * (roh**2)
    res = ne/za
    tau = Q * (res)**(1./3.)
    return tau

def transmission(tau, theta):
    """
    Diese Funktion berechnet die Transmission bei gegebener optischen Tiefe
    und Zenitwinkel
    
    Input:
        tau : optische Dicke in [/]
        theta : zenitwinkel in [°]
        
    Output:
        transmission in [/]
    """
    
    res = -1. * (tau / np.cos(np.deg2rad(theta)))
    t = np.exp(res)
    return t
              
    
def swimmingpool(lam, ni, x):
    """
    Funktion zur Berechnung der Transmission im Wasser bei gegebener
    Wellenlaenge und imaginaeren Brechungsindex.
    
    Input:
        lam : Waellenlaenge in mikrometer
        ni : imaginaerteil Brechungsindex dimensionslos
        x : Abstand in mikrometern
    Output:
        t : Transmission dimensionslos
    """
    ne = -1 * (4 * np.pi * ni * x) 
    ze = lam
    res = ne / ze
    t = np.exp(res)
    return t


def Tb_sat(lam, intens):
    """
    Diese Funktion berechnet aus gegebener Wellenlaenge und
    Strahlungsintensitaet die Strahlungstemperatur über die invetierte
    Planck Funktion.
    
    Input:
        I : Intensitaet in [W m^-2 µm^-1 sr^-1 ]
            (mal 10e6 um die Micrometer verrechnen zu können)
        lam : Wellenlaenge in [µm]
            (Wellenlaenge in Micrometer)
        
    Output:
        tb : Strahlungstemperatur 
    
    """
    
    a1 = (2.0 * h * c**2.) / (intens * lam**5.) 
    #print (a1)
    b1 = (h * c) / (k * lam)
    #print (b1)
    tb = b1 * 1./np.log(1.0 + a1)
    return tb


def size_para(rr,ww):
    """
    Diese Funktion berechnet den Size Parameter in Abhaengigkeit von
    Partikelradius zu Wellenlaenge. GP (12.1)
    
    Input: 
        rr : Partikelradius in Micrometer [µm]
        ww : Wellenlaenge in Micrometer [µn]
    Output:
        chi : Size Parameter ohne Dimension [/]
    """
    chi = (2. * np.pi * rr) / ww
    return chi

def ray_efficiencies(m,chi):
    """
    Diese Funktion berechnet aus dem Brechungsindex und dem Size Parameter die
    Efficiensies mit der Rayleigh Theorie. (GP 12.12 - 12.14)

    Input:
        m : Brechungsindex in [/]
        chi : Size Parameter in [/]
        
    Output:
        Qe : extinktions Efficiency
        Qs : streuungs Efficiency
        Qa : absorbtions Efficiency
    """
    res = 1. + (chi*+2./15.)*((m**2. -1.) / (m**2. + 2.))*((m**4. + 27.*m**2.+38.)/(2.*m**2. +3.))
    Qe = 4.*chi*(((m**2. -1.) / (m**2. + 2.))*res).imag
    Qs = (8./3.) * chi**4. * abs((m**2. -1.) / (m**2. + 2.))**2.
    Qa =4.*chi*((m**2. -1.) / (m**2. + 2.)).imag 
    
    return Qe, Qs, Qa


def efficiencies(m,chi):
    """
    Diese Funktion berechnet aus dem Brechungsindex und dem Size Parameter die
    Efficiensies mit der Mie Theorie.


    Input:
        m : Brechungsindex in [/]
        chi : Size Parameter in [/]

    Output:
        Qe : extinktions Efficiency
        Qs : streuungs Efficiency
        Qa : absorbtions Efficiency
    """
    
    ## calc N
    N=round(chi+4.*chi**(1./3.)+2,0)
    #N=chi+4*chi**(1/3)+2
    
    ## initialize Mie coefficients
    an=np.zeros(int(N), dtype=complex)
    bn=np.zeros(int(N), dtype=complex)
    
    ## initialize efficiencies
    Qext, Qsca, Qabs = 0, 0, 0

    ## initialize recursion terms for Mie coefficients
    W0=complex(np.sin(chi), np.cos(chi))
    Wm1=complex(np.cos(chi),-np.sin(chi))
    A0=1./(np.tan(m*chi))
    ## loop over Mie coefficients and sum contributions to
    ## the efficiencies
    n=0
    while(n<N):
        n = n+1
        An = -n/(m*chi)+1/(n/(m*chi)-A0)
        Wn = ((2*n-1)/chi)*W0-Wm1
        an[n-1] = ((An/m+n/chi)*Wn.real-W0.real)/((An/m+n/chi)*Wn-W0)
        bn[n-1] = ((m*An+n/chi)*Wn.real-W0.real)/((m*An+n/chi)*Wn-W0)
        A0 = An
        Wm1 = W0
        W0 = Wn
        Qext = Qext+2/(chi**2)*(2*n+1)*np.real(an[n-1]+bn[n-1])
        Qsca = Qsca+2/(chi**2)*(2*n+1)*(abs(an[n-1])**2+abs(bn[n-1])**2)

    ## globalize Mie coefficients and efficiencies
    return Qext, Qsca, Qext-Qsca, an, bn


def ray_phase_func(theta):
    """
    Diese Funktion berechnet aus dem Streuwinkel die Pahsenfunktion nach
    der Rayleigh Theorie. GP (12.10)

    Input:
        theta : Streuwinkel in [°]

    Output:
        p : Rayleigh Phasen Funktion in [/]
    """
    p = (3./4.) *(1.+np.cos(np.deg2rad(theta))**2)
    return p



def phase_func(m,chi,mu, nang):
    """
    Diese Funktion berechnet aus dem Streuwinkel die Pahsenfunktion nach
    der Mie Theorie. ()

    Input:
        m : Brechungsindex in [/]
        chi : Size Parameter in [/]
        mu : Cosinus des Streuwinkels in [/]
        nang : Länge des Streuwinkelvektors in [/]

    Output:
        p : Mie Phasen Funktion in [/]
    """
    
    N = np.round(chi + 4.*chi**(1./3.) + 2.,0)
    
    Qext, Qsca, Qabs, an, bn = efficiencies(m,chi)

    factor=2./(chi*chi*Qsca)

    S1 = np.zeros(nang, dtype=complex)
    S2 = np.zeros(nang, dtype=complex)
    p11 = np.zeros(nang, dtype=complex)

    i=0

    while(i<(nang)):

        i=i+1
        n=1
        pi0 = 0
        pi1 = 1
        tau1 = mu[i-1]

        while(n<N):
            S1[i-1]=S1[i-1]+(2*n+1)/(n*(n+1))*(an[n-1]*pi1+bn[n-1]*tau1)
            S2[i-1]=S2[i-1]+(2*n+1)/(n*(n+1))*(bn[n-1]*pi1+an[n-1]*tau1)
            n=n+1
            pi1n=(2*n-1)/(n-1)*mu[i-1]*pi1-n/(n-1)*pi0
            tau1=n*mu[i-1]*pi1n-(n+1)*pi1
            pi0=pi1
            pi1=pi1n

        p11[i-1]=factor*(S1[i-1]*np.conj(S1[i-1])+S2[i-1]*np.conj(S2[i-1]))
    
    return p11



def weg(a,rando):
    """
    Diese Funktion berechnet aus dem Volumenextinktionskoeffizienten und einer
    Zufallszahl die Weglänge die ein Photon in einem Medium passiert.
    

    Input:
        a : Volumenextinktionskoeffizient in [1/m]
        rando : Zufallszahl in [/]
        
    Output:
        l : Weglänge in [m]
        
    """    
    l = -1. * np.log(1-rando) / a # 
    return l


def Theta_HG(random,g):
    """
    Diese Funktion brechnet den zufälligen Streuwinkel nach Henyey Greenstein
    unter Angabe des Asymetrieparameters und einer Zufallszahl. (GP 11.23)
        
    Input:
        random : Zufallszahl in [/]
        g : Asymetrieparameter in [/]

    Output:
        Theta : Streuwinkel in [rad]
        
    """    
    if g==0:
        Theta=np.arccos(2*random-1);
        return Theta
    else:
        #par=(2. * g * random / (1. - g**2) + 1. / (1. +g ))**2;
        #mu=(1+g**2-1./par)/(2*g);
        #Theta = np.arccos(mu)
        res = (1. + g**2) - ((1.- g**2.)/(1.- g + 2.*g*random))**2
        Theta = np.arccos((1./(2. *g)*res))
        return Theta


def rotmat(theta_old,phi_old):
    """
    Diese Funktion berechnet die Rotationmatrix der Richtung aus dem
    Einheitssystem in das Absolutsystem. ()

    Input:
        theta_old : vorheriger Elevationswinkel
        phi_old : vorheriger Azimuthwinkel
        
    Output:
        3x3 Rotationsmatrix
        
    """    
    R11 = np.cos(theta_old) * np.cos(phi_old)
    R21 = np.cos(theta_old) * np.sin(phi_old)
    R31 = -np.sin(theta_old)
    R12 = -np.sin(phi_old)
    R22 = np.cos(phi_old)
    R32 = 0
    R13 = np.sin(theta_old) * np.cos(phi_old)
    R23 = np.sin(theta_old) * np.sin(phi_old)
    R33 = np.cos(theta_old)
    A = np.matrix([[R11,R21,R31], [R12,R22,R32], [R13,R23,R33]])
    return np.matrix.transpose(A)
    


def H_von_OptDic(beta_e, opt_dic):
    """
    Berechnung der Höhe bei gegebener optischen Dicke und einem
    Volumenextinktionskoeffizient.
    GP (7.32)
        
    Input:
        beta_e : Volumenextinktionskoeffizient in [1/m]
        opt_dic: optische Dicke in [/]
        
    Output:
        H : Höhe der Schicht in [m]
        
    """    

    
    H = opt_dic/beta_e
    
    return H


def gamma(omega_bar, g):
    """
    Berechnung von Gamma als Zwischenergebnis der Zweistromapproximation
    GP (13.25)

    Input:
        omega_bar : Einfachstreualbedo in [/]
        g: Aysmetrieparameter in [/]

    Output:
        gamma : Gamma in [/]

    """
    gamma = 2*np.sqrt(1-omega_bar)*np.sqrt(1-omega_bar*g)
    return gamma

def r_inf(omega_bar, g):
    """
    Berechnung der Albedo am Oberand der Wolke und uter der Annahme der
    Zweistromapproximation bei einer semi-unendlichen Wolke.
    GP (13.44)

    Input:
        omega_bar : Einfachstreualbedo in [/]
        g: Aysmetrieparameter in [/]

    Output:
        r_inf

    """
    r_inf = (np.sqrt(1. - omega_bar * g) - np.sqrt(1. - omega_bar)) /\
            (np.sqrt(1. - omega_bar * g) + np.sqrt(1. - omega_bar))
    return r_inf


def tsa_r(g,tau_stern, omega_bar):
    """
    Berechnung der Reflexion am Oberand der Wolke und unter der Annahme der
    Zweistromapproximation .
    GP (13.65)

    Input:
        omega_bar : Einfachstreualbedo in [/]
        g: Aysmetrieparameter in [/]
        tau_stern : optische Dicke in [/]

    Output:
        r : Reflexion in [/]

    """
    if omega_bar == 1:
        r = ((1. - g) * tau_stern) / (1. + (1. - g) * tau_stern)
    else:
        n = r_inf(omega_bar, g) * \
            (np.exp(gamma(omega_bar, g) * tau_stern) -
             np.exp(-1. * gamma(omega_bar, g) * tau_stern))

        z = np.exp(gamma(omega_bar, g) *
                   tau_stern) - np.exp(-1. * gamma(omega_bar, g) * tau_stern)\
                                * r_inf(omega_bar, g)**2
        r = n / z
    return r

def tsa_t(g,tau_stern, omega_bar):
    """
    Berechnung der Transmission am Oberand der Wolke und unter der Annahme der
    Zweistromapproximation .
    GP (13.66)

    Input:
        omega_bar : Einfachstreualbedo in [/]
        g: Aysmetrieparameter in [/]
        tau_stern : optische Dicke in [/]

    Output:
        t : Reflexion in [/]

    """
    if omega_bar == 1:
        
        t = (1.) / (1. + (1. - g) * tau_stern)
        
    else:
        n = 1 - r_inf(omega_bar, g)**2
        z = np.exp(gamma(omega_bar, g) * tau_stern) - \
                np.exp(-1. * gamma(omega_bar, g)*tau_stern) * \
                r_inf(omega_bar, g)**2
        t = n / z
    return t

def tsa_tdiff(g,tau_stern, omega_bar):
    """
    Berechnung der direkten und diffusen Transmission einer Wolke unter der
    Annahme der Zweistromapproximation .
    GP (13.69)

    Input:
        omega_bar : Einfachstreualbedo in [/]
        g: Aysmetrieparameter in [/]
        tau_stern : optische Dicke in [/]

    Output:
        tdiff : Diffuse transmission in [/]

    """
    if omega_bar == 0:
        tdiff = 0
    if omega_bar == 1:
        tdiff = ((1.) / (1. + (1. - g) * tau_stern)) - np.exp(-tau_stern / 0.5)
    else:
        n = 1-r_inf(omega_bar, g)**2
        z = z = np.exp(gamma(omega_bar, g) * tau_stern) - \
                np.exp(-1. * gamma(omega_bar, g) * tau_stern) *\
                r_inf(omega_bar, g)**2

        tdiff = (n / z) - np.exp(-tau_stern / 0.5)
        
    return tdiff
