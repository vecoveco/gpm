#!/usr/bin/env python

#WEITERES MODUL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#https://github.com/jleinonen/pytmatrix/wiki#scatter-module


# Author: Benjamin Marchant
# Ref: Hong Du's paper: "Mie-Scattering Calculation, Applied Optics, Vol. 43, Issue 9, pp. 1951-1956 (2004)

import matplotlib.pyplot as plt

import numpy as np
import math

#----------------------------------------------------------------------------------------#
# user inputs

ri = complex(1.33,-0.0) # refractive index

x = 5.0 # Mie Parameter (2*pi*radius/wavelength)

n = 100

angular_resolution = 1.0

#----------------------------------------------------------------------------------------#
# Calculations of the Riccati-Bessel functions by upward recurrence.

PHI = np.zeros((n+1))

a = math.cos(x)
PHI[0] = math.sin(x)
b = PHI[0]
for i in range(1,n+1):
	PHI[i] = ( ( 2 * i - 1 ) / x ) * b  - a
	a = b
	b = PHI[i]

#----------------------------------------------------------------------------------------#

CHI = np.zeros((n+1))

a = math.cos(x)
CHI[0] = a
CHI[1] = math.cos(x) / x + math.sin(x)
b = CHI[1]
for i in range(2,n+1):
	CHI[i] = ( ( 2 * i - 1 ) / x ) * b  - a
	a = b
	b = CHI[i]

#----------------------------------------------------------------------------------------#

ZETA = np.zeros((n+1), dtype=complex)
for i in range(0,n+1):
	ZETA[i] = complex(PHI[i], CHI[i])

#----------------------------------------------------------------------------------------#
# Calculations of the complex functions rn(mx) by upward recurrence.

RN = np.zeros((n+1), dtype=complex)

a = ri.real * x
b = ri.imag * x

c1 = math.tan(a) - math.exp(-2.0*b) * math.tan(a)
c2 = 1.0 + math.exp(-2.0*b)
c3 = -1.0 + math.exp(-2.0*b)
c4 = math.tan(a) + math.exp(-2.0*b) * math.tan(a)

a = c1 * c3 + c2 * c4
b = c2 * c3 - c1 * c4
c = c3 * c3 + c4 * c4

RN[0] = complex(a/c, b/c)
for i in range(1,n+1):
	a = ( 2.0 * i - 1.0 ) / ( ri * x ) - RN[i-1]
	RN[i] = 1.0 / a

#----------------------------------------------------------------------------------------#
# Calculations of the Mie coefficients an et bn.

MA = np.zeros((n), dtype=complex)
MB = np.zeros((n), dtype=complex)

for i in range(0,n):
	a = RN[i+1] / ri + (i + 1) * ( 1.0 - 1.0 / ( ri * ri ) ) / x
	MA[i] = ( a * PHI[i+1] - PHI[i]) / ( a * ZETA[i+1] - ZETA[i] )
	MB[i] = ( RN[i+1] * ri * PHI[i+1] - PHI[i] ) / ( RN[i+1] * ri * ZETA[i+1] - ZETA[i] )

#----------------------------------------------------------------------------------------#

def angular_functions(x,angle,n,MA,MB):

	#------------------------------------------------------------------------------------#
	# Calculations of the Mie angular functions.

	PI = np.zeros((n))

	a = math.cos(angle)
	b = 0.0
	c = 1.0

	PI[0] = 1.0
	for j in range(1,n):
		s = a * c
		t = s - b
		PI[j] = s + t + t / j
		b = c
		c = PI[j]

	#------------------------------------------------------------------------------------#

	TAU = np.zeros((n))

	TAU[0] = a
	for j in range(1,n):
		t = a * PI[j] - PI[j-1]
		TAU[j] = ( j + 1 ) * t - PI[j-1]

	#------------------------------------------------------------------------------------#
	# scattering amplitudes

	b = complex(0.0,0.0)
	c = complex(0.0,0.0)
	for i in range(0,n):
		a = ( 2.0 * ( i + 1 ) + 1.0 ) / ( ( i + 1 ) * ( i + 2 ) )
		b = b + a * ( MA[i] * PI[i] + MB[i] * TAU[i] )
		c = c + a * ( MB[i] * PI[i] + MA[i] * TAU[i] )
	sa1 = b
	sa2 = c

	return sa1, sa2

#----------------------------------------------------------------------------------------#
# phase matrix

angle_list = []
p11_list = []
p12_list = []
p33_list = []
p34_list = []

for angle in np.arange(0.0,180.0,angular_resolution):

	#print angle

	sa1,sa2 = angular_functions(x,math.radians(angle),n,MA,MB)

	p11 = ( sa1 * sa1.conjugate() + sa2 * sa2.conjugate() )
	p12 = ( sa2 * sa2.conjugate() - sa1 * sa1.conjugate() )
	p33 = ( sa2 * sa1.conjugate() + sa2 * sa1.conjugate() )
	p34 = complex(0,1) * ( sa1 * sa2.conjugate() - sa2 * sa1.conjugate() )

	angle_list.append(angle)
	p11_list.append(p11.real)
	p12_list.append(p12.real)
	p33_list.append(p33.real)
	p34_list.append(p34.real)

#----------------------------------------------------------------------------------------#
# normalization

norm = 0

for i in range(0,n):
	a = abs(MA[i]) * abs(MA[i]) + abs(MB[i]) * abs(MB[i])
	norm = norm +  ( 2.0 * ( i + 1.0 ) + 1.0 ) * a

p11_list = p11_list / norm
p12_list = p12_list / norm
p33_list = p33_list / norm
p34_list = p34_list / norm

#----------------------------------------------------------------------------------------#
# extinction efficiency.

a = 0.0
for i in range(0,n):
	a = a + ( 2 * ( i + 1 ) + 1 ) * ( MA[i] + MB[i] ).real

qe = ( 2.0 / ( x * x ) ) * a

#----------------------------------------------------------------------------------------#
# asymmetry factor.

a = 0.0
for i in range(0,n-1):
	c1 = MA[i] * MA[i+1].conjugate() + MB[i] * MB[i+1].conjugate()
	c2 = MA[i] * MB[i].conjugate()
	b = ( (i+1) * ( (i+1) + 2.0 ) ) / ( (i+1) + 1.0 )
	c = ( 2.0 * (i+1) + 1.0 ) / ( (i+1) * ( (i+1) + 1.0 ) )
	a = a + b * c1.real + c * c2.real

gfact = ( 2.0 / norm ) * a

#----------------------------------------------------------------------------------------#
# single scattering albedo.

albedo = ( 2.0 / ( x * x ) ) * norm / qe

#----------------------------------------------------------------------------------------#
# Data Visualization

fig = plt.figure()
ax = fig.add_subplot(111)

plt.plot(angle_list,p11_list,'k-')
plt.yscale('log')

plt.grid(True,which="both")

plt.xlabel(r"Scattering Angle $\Theta$ ($^\circ$)")
plt.ylabel(r"$P_{11}$")

plt.text(0.70,0.85,'x = '+str(x), fontsize=12, transform = ax.transAxes)
plt.text(0.70,0.80,'n=('+str(ri.real)+','+str(ri.imag)+')', \
                    fontsize=12, transform = ax.transAxes)

plt.text(0.70,0.70,'qe = '+str(round(qe,2)), fontsize=12, transform = ax.transAxes)
plt.text(0.70,0.65,'w0 = '+str(round(albedo,2)), fontsize=12, transform = ax.transAxes)
plt.text(0.70,0.60,'g = '+str(round(gfact,2)), fontsize=12, transform = ax.transAxes)

#plt.savefig('p11.png')
plt.show()

#----------------------------------------------------------------------------------------#
# plot p12

fig = plt.figure()
ax = fig.add_subplot(111)

plt.plot(angle_list,p12_list,'k-')
#plt.yscale('log')

plt.grid(True,which="both")

plt.xlabel(r"Scattering Angle $\Theta$ ($^\circ$)")
plt.ylabel(r"$P_{12}$")

#plt.savefig('p12.png')
plt.show()

#----------------------------------------------------------------------------------------#
# plot p33

fig = plt.figure()
ax = fig.add_subplot(111)

plt.plot(angle_list,p33_list,'k-')
plt.yscale('log')

plt.grid(True,which="both")

plt.xlabel(r"Scattering Angle $\Theta$ ($^\circ$)")
plt.ylabel(r"$P_{33}$")

#plt.savefig('p33.png')
plt.show()

#----------------------------------------------------------------------------------------#
# plot p34

fig = plt.figure()

ax = fig.add_subplot(111)

plt.plot(angle_list,p34_list,'k-')
#plt.yscale('log')

plt.grid(True,which="both")

plt.xlabel(r"Scattering Angle $\Theta$ ($^\circ$)")
plt.ylabel(r"$P_{34}$")

#plt.savefig('p34.png')
plt.show()






################################################################

"""
Copyright (C) 2012-2013 Jussi Leinonen
Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""


from numpy import arange, dot, zeros, vstack, sqrt, sin, cos


def mie_props(coeffs,y):
    """The scattering properties.
    """
    anp = coeffs.an.real
    anpp = coeffs.an.imag
    bnp = coeffs.bn.real
    bnpp = coeffs.bn.imag
    nmax = coeffs.nmax

    n1 = nmax-1
    n = arange(1,nmax+1,dtype=float)
    cn = 2*n+1
    c1n = n*(n+2)/(n+1)
    c2n = cn/(n*(n+1))
    y2 = y**2

    dn = cn*(anp+bnp)
    q = dn.sum()
    qext = 2*q/y2

    en = cn*(anp**2+anpp**2+bnp**2+bnpp**2)
    q = en.sum()
    qsca = 2*q/y2
    qabs = qext-qsca

    fn = (coeffs.an-coeffs.bn)*cn
    gn=(-1)**n
    q = (fn*gn).sum()
    qb = dot(q,q.conj()).real/y2

    g1 = zeros((4,nmax),dtype=float)
    g1[:,:n1] = vstack((anp[1:nmax], anpp[1:nmax], bnp[1:nmax], bnpp[1:nmax]))

    asy1 = c1n*(anp*g1[0,:]+anpp*g1[1,:]+bnp*g1[2,:]+bnpp*g1[3,:])
    asy2 = c2n*(anp*bnp+anpp*bnpp)

    asy = 4/y2*(asy1+asy2).sum()/qsca
    qratio = qb/qsca

    return {"qext":qext, "qsca":qsca, "qabs":qabs, "qb":qb, "asy":asy,
        "qratio":qratio}


def mie_S12(coeffs,u):
    """The amplitude scattering matrix.
    """
    (pin,tin) = mie_pt(u,coeffs.nmax)
    n = arange(1, coeffs.nmax+1, dtype=float)
    n2 = (2*n+1)/(n*(n+1))
    pin *= n2
    tin *= n2

    S1 = dot(coeffs.an,pin)+dot(coeffs.bn,tin)
    S2 = dot(coeffs.an,tin)+dot(coeffs.bn,pin)
    return (S1, S2)


def mie_pt(u,nmax):
    u = float(u)
    p = zeros(nmax, dtype=float)
    p[0] = 1
    p[1] = 3*u
    t = zeros(nmax, dtype=float)
    t[0] = u
    t[1] = 6*u**2 - 3

    nn = arange(2,nmax,dtype=float)
    for n in nn:
        n_i = int(n)
        p[n_i] = (2*n+1)/n*p[n_i-1]*u - (n+1)/n*p[n_i-2]

    t[2:] = (nn+1)*u*p[2:] - (nn+2)*p[1:-1]


	return (p,t)


"""
Copyright (C) 2012-2013 Jussi Leinonen
Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from numpy import pi, arange, zeros, hstack, sqrt, sin, cos
from scipy.special import jv, yv


class MieCoeffs(object):
    """Wrapper for the Mie coefficients.
    """
    def __init__(self, par):
        (self.an, self.bn, self.nmax) = mie_coeffs(par)


def mie_coeffs(params):
    """Input validation and function selection for the Mie coefficients.
    """

    eps = complex(params["eps"]) if params["eps"] is not None else None
    x = float(params["x"]) if params["x"] is not None else None

    if (x==None) or (eps==None):
        raise ValueError("Must specify x and either eps or m.")

    mu = complex(params["mu"]) if params["mu"] is not None else complex(1.0)

    y = float(params["y"]) if params["y"] is not None else None
    eps2 = complex(params["eps2"]) if params["eps2"] is not None else None

    coated = (y is not None)
    if coated == (eps2 is None):
        raise ValueError("Must specify both y and m2 for coated particles.")

    if coated and mu!=complex(1.0):
        raise ValueError("Multilayer calculations for magnetic particles " + \
            "are not currently supported.")

    if not coated:
        y = x
        eps2 = eps

    # Do not use the coated version if it is not necessary
    if x==y or eps==eps2:
        coeffs = single_mie_coeff(eps,mu,y)
    elif x==0:
        coeffs = single_mie_coeff(eps2,mu,y)
    else:
        coeffs = coated_mie_coeff(eps,eps2,x,y)

    return coeffs


def single_mie_coeff(eps,mu,x):
    """Mie coefficients for the single-layered sphere.
    Args:
        eps: The complex relative permittivity.
        mu: The complex relative permeability.
        x: The size parameter.
    Returns:
        A tuple containing (an, bn, nmax) where an and bn are the Mie
        coefficients and nmax is the maximum number of coefficients.
    """
    z = sqrt(eps*mu)*x
    m = sqrt(eps/mu)

    nmax = int(round(2+x+4*x**(1.0/3.0)))
    nmax1 = nmax-1
    nmx = int(round(max(nmax,abs(z))+16))
    n = arange(nmax)
    nu = n+1.5

    sx = sqrt(0.5*pi*x)
    px = sx*jv(nu,x)
    p1x = hstack((sin(x), px[:nmax1]))
    chx = -sx*yv(nu,x)
    ch1x = hstack((cos(x), chx[:nmax1]))
    gsx = px-complex(0,1)*chx
    gs1x = p1x-complex(0,1)*ch1x

    dnx = zeros(nmx,dtype=complex)
    for j in range(nmx-1,0,-1):
        r = (j+1.0)/z
        dnx[j-1] = r - 1.0/(dnx[j]+r)
    dn = dnx[:nmax]
    n1 = n+1
    da = dn/m + n1/x
    db = dn*m + n1/x

    an = (da*px-p1x)/(da*gsx-gs1x)
    bn = (db*px-p1x)/(db*gsx-gs1x)

    return (an, bn, nmax)


def coated_mie_coeff(eps1,eps2,x,y):
    """Mie coefficients for the dual-layered (coated) sphere.
       Args:
          eps: The complex relative permittivity of the core.
          eps2: The complex relative permittivity of the shell.
          x: The size parameter of the core.
          y: The size parameter of the shell.
       Returns:
          A tuple containing (an, bn, nmax) where an and bn are the Mie
          coefficients and nmax is the maximum number of coefficients.
    """
    m1 = sqrt(eps1)
    m2 = sqrt(eps2)
    m = m2/m1
    u = m1*x
    v = m2*x
    w = m2*y

    nmax = int(round(2+y+4*y**(1.0/3.0)))
    mx = max(abs(m1*y),abs(w))
    nmx = int(round(max(nmax,mx)+16))
    nmax1 = nmax-1
    n = arange(nmax)

    dnu = zeros(nmax,dtype=complex)
    dnv = zeros(nmax,dtype=complex)
    dnw = zeros(nmax,dtype=complex)
    dnx = zeros(nmx,dtype=complex)

    for (z, dn) in zip((u,v,w),(dnu,dnv,dnw)):
        for j in range(nmx-1,0,-1):
            r = (j+1.0)/z
            dnx[j-1] = r - 1.0/(dnx[j]+r)
        dn[:] = dnx[:nmax]

    nu = n+1.5
    vwy = [v,w,y]
    sx = [sqrt(0.5*pi*xx) for xx in vwy]
    (pv,pw,py) = [s*jv(nu,xx) for (s,xx) in zip(sx,vwy)]
    (chv,chw,chy) = [-s*yv(nu,xx) for (s,xx) in zip(sx,vwy)]
    p1y = hstack((sin(y), py[:nmax1]))
    ch1y = hstack((cos(y), chy[:nmax1]))
    gsy = py-complex(0,1)*chy
    gs1y = p1y-complex(0,1)*ch1y

    uu = m*dnu-dnv
    vv = dnu/m-dnv
    fv = pv/chv
    fw = pw/chw
    ku1 = uu*fv/pw
    kv1 = vv*fv/pw
    pt = pw-chw*fv
    prat = pw/pv/chv
    ku2 = uu*pt+prat
    kv2 = vv*pt+prat
    dns1 = ku1/ku2
    gns1 = kv1/kv2

    dns = dns1+dnw
    gns = gns1+dnw
    nrat = (n+1)/y
    a1 = dns/m2+nrat
    b1 = m2*gns+nrat
    an = (py*a1-p1y)/(gsy*a1-gs1y)
    bn = (py*b1-p1y)/(gsy*b1-gs1y)

	return (an, bn, nmax)