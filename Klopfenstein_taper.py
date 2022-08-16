"""
@author: Ziad (https://github.com/ZiadHatab)
@email: zi.hatab@gmail.com
under BSD 3-Clause License
"""
import numpy as np  # python -m pip install -U numpy
import functools    # standard library
import matplotlib.pyplot as plt  # python -m pip install -U matplotlib

'''
# References
[1] R. W. Klopfenstein, 
"A Transmission Line Taper of Improved Design," 
in Proceedings of the IRE, vol. 44, no. 1, pp. 31-35, Jan. 1956, doi: 10.1109/JRPROC.1956.274847.
https://ieeexplore.ieee.org/document/4051841

[2] D. Kajfez and J. O. Prewitt, 
"Correction to "A Transmission Line Taper of Improved Design" (Letters)," 
in IEEE Transactions on Microwave Theory and Techniques, vol. 21, no. 5, pp. 364-364, May 1973, doi: 10.1109/TMTT.1973.1128003.
https://ieeexplore.ieee.org/document/1128003

[3] M. A. Grossberg, 
"Extremely rapid computation of the Klopfenstein impedance taper," 
in Proceedings of the IEEE, vol. 56, no. 9, pp. 1629-1630, Sept. 1968, doi: 10.1109/PROC.1968.6686.
https://ieeexplore.ieee.org/document/1448616

[4] Michael Steer, 
Microwave and RF Design: Networks. Volume 3. (Third Edition), 
NC StateUniversity, 2019. doi: https://doi.org/10.5149/9781469656953_Steer
'''

def phi(z,A=1,N=40):
    '''
    The phi function defined by Klopfenstein in [1].
    Implementation based on finite sum approximation using the procedure in [3].
    '''
    an = 1
    bn = z/2
    pn = [an*bn]
    for n in np.arange(N-1)+1:
        an = A**2*an/4/n/(n+1)
        bn = (z/2*(1-z**2)**n + 2*n*bn)/(2*n + 1)
        pn.append(an*bn)
    return np.sum(pn)
'''
import scipy.special as ss
import scipy.integrate as si
def phi(z,A=1):
    # directly compute the integral with scipy integral method and Bessel function.
    return si.quad(lambda y: ss.iv(1, A*np.sqrt(1-y**2))/A/np.sqrt(1-y**2), 0, z)[0]
'''

def klopf(Z1, Z2, Gmax, N=100):
    '''
    Klopfenstein transmission line taper [1].
    Check [4] for a text book explanation on Klopfenstein taper.
    Z1  : start impedance (scalar)
    Z2  : end impedance (scalar)
    Gmax: maximum reflection in passband (scalar)
    N   : number of segments points (scalar)
    '''
    U = lambda x: np.heaviside(x, 1)  # step function
    lnZ1Z2 = np.log(Z1*Z2)
    G0     = np.log(Z2/Z1)/2
    A      = np.arccosh(abs(G0/Gmax))
    L      = 1  # length of the taper. It can be anything > 0, because all units are normalized in below equation.
    # I modified the correction of [2] so that Z starts with Z1 and ends with Z2 without needing to split the equation.
    lnZ = np.array([lnZ1Z2/2 + G0/np.cosh(A)*(A**2*phi(2*x/L,A) + U(x-L/2) - U(-x-L/2)) for x in np.linspace(-L/2, L/2, N)])
    return np.exp(lnZ)

def klopf_S(Z1, Z2, Gmax, L, f, ereff=1-0j):
    '''
    The theoretical S-parameters of the Klopfenstein taper (Tchebycheff filter).
    Z1   : start impedance (scalar)
    Z2   : end impedance (scalar)
    Gmax : maximum reflection in passband (scalar)
    L    : length of the taper in meters (scalar)
    f    : frequency points in Hz (array)
    ereff: relative effective permittivity (scalar)
    '''
    c0  = 299792458   # speed of light in vacuum (m/s)
    g   = 2*np.pi*f*np.sqrt(-complex(ereff))/c0  # propagation constant
    g   = g*np.sign(g)        # correct for the sign convention of ereff = real - 1j*imag 
    G0  = (Z2 -Z1)/(Z2 + Z1)  # impedance mismatch 
    A   = np.arccosh(abs(G0/Gmax))
    S11 = np.exp(-g*L)*G0*np.cosh(np.sqrt((g*L)**2 + A**2))/np.cosh(A) # modified for lossy response
    S21 = np.sqrt(1-abs(S11)**2)*np.exp(-g*L)
    return np.array([ [[s11,s21],[s21,-s11]] for s11,s21 in zip(S11,S21) ])

def klopf_f2L(Z1, Z2, Gmax, f3db, ereff=1):
    '''
    Determining the minimum length of the Klopfenstein taper that correspond to a given 3dB cut-off frequency.
    Z1   : start impedance (scalar)
    Z2   : end impedance (scalar)
    Gmax : maximum reflection in passband (scalar)
    f3db : 3dB cut-off frequency in Hz (scalar)
    ereff: relative effective permittivity(scalar)
    '''
    c0  = 299792458   # speed of light in vacuum (m/s)
    G0  = (Z2 -Z1)/(Z2 + Z1)  # impedance mismatch 
    A   = np.arccosh(abs(G0/Gmax))
    return c0*np.sqrt(np.arccos(1/2)**2 + A**2)/2/np.pi/np.sqrt(ereff.real)/f3db

def klopf_L2f(Z1, Z2, Gmax, L, ereff=1):
    '''
    Determining the low 3dB cut-off frequency of Klopfenstein taper given its length.
    Z1   : start impedance (scalar)
    Z2   : end impedance (scalar)
    Gmax : maximum reflection in passband (scalar)
    L    : length of the taper in meters (scalar)
    ereff: relative effective permittivity (scalar)
    '''
    c0  = 299792458   # speed of light in vacuum (m/s)
    G0  = (Z2 -Z1)/(Z2 + Z1)  # impedance mismatch 
    A   = np.arccosh(abs(G0/Gmax))
    return c0*np.sqrt(np.arccos(1/2)**2 + A**2)/2/np.pi/np.sqrt(ereff.real)/L

def S2T(S):
    # convert S to T parameters
    T = S.copy()
    T[0,0] = -(S[0,0]*S[1,1]-S[0,1]*S[1,0])
    T[0,1] = S[0,0]
    T[1,0] = -S[1,1]
    T[1,1] = 1
    return T/S[1,0]

def T2S(T):
    # convert T to S parameters
    S = T.copy()
    S[0,0] = T[0,1]
    S[0,1] = T[0,0]*T[1,1]-T[0,1]*T[1,0]
    S[1,0] = 1
    S[1,1] = -T[1,0]
    return S/T[1,1]

def Qnm(Zn, Zm):
    '''
    Impedance transformer in T-parameters from Zn to Zm. Based on the Eqs. (86) and (87) in
    R. Marks and D. Williams, "A general waveguide circuit theory," 
    Journal of Research (NIST JRES), National Institute of Standards and Technology,
    Gaithersburg, MD, no. 97, 1992.
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4914227/
    
    Zn : start impedance (scalar)
    Zm : end impedance (scalar)
    '''
    Gnm = (Zm-Zn)/(Zm+Zn)
    return np.sqrt(Zn.real/Zm.real*(Zm/Zn).conjugate())/np.sqrt(1-Gnm**2)*np.array([[1, Gnm],[Gnm, 1]])

def taperS(Z, L, f, ereff=1-0j):
    '''
    S-parameters of an arbitrary impedance taper defined by its Z profile and length.
    Implemented as cascaded segments using T-parameters. Every two values in Z is one segment. 
    Z    : impedance profile (array)
    L    : length of the taper in meters (scalar)
    f    : frequency points in Hz (array)
    ereff: relative effective permittivity (scalar)
    '''
    dT = lambda x,Zn,Zm: Qnm(Zn, Zm)@np.diag([np.exp(-x), np.exp(x)])  # T-parameters of one segment
    c0 = 299792458   # speed of light in vacuum (m/s)
    g  = 2*np.pi*f*np.sqrt(-complex(ereff))/c0  # propagation constant
    g  = g*np.sign(g)   # correct for the sign convention of ereff = real - 1j*imag 
    dl = L/(len(Z)-1)   # length of a single segment    
    return np.array([T2S( functools.reduce(np.dot, [dT(p*dl,z1,z2) for z1,z2 in zip(Z[:-1],Z[1:])]) ) for p in g])
    

if __name__ == '__main__':
    mag2db = lambda x: 20*np.log10(abs(x))
    db2mag = lambda x: 10**(x/20)
    
    Z1 = 50
    Z2 = 120
    N  = 200     # number of segments 
    # the more segments, the better. However, it gets slower when computing S-parameters by matrix cascade.
    Gmax  = db2mag(-40) # max reflection in passband
    # Gmax must be lower than (Z2-Z1)/(Z2+Z1), otherwise you are introducing gain in S11. I will let you think about that!
    Zklopf  = klopf(Z1, Z2, Gmax, N)  # Klopfenstein impedance profile
    
    x = np.linspace(0,1,N) # normlized length
    # plot the impedance profile
    fig, ax = plt.subplots(1,1, figsize=(5.5, 5.5/1.5))
    ax.plot(x, Zklopf, '-' , lw=2, label='Klopfenstein')
    ax.set_ylabel('Impedance profile (ohm)')
    ax.set_xlabel('Normalized length')
    fig.set_dpi(600)
    fig.tight_layout(pad=1.2)
    plt.legend()
    
    # plot the frequency response (theoretical vs finite cascaded segments)
    M = 400   # number of frequency points
    f = np.linspace(1,150,M)*1e9 # in Hz
    L = 4e-3 # in meter
    ereff = 2.5 - 0j
    Sklopf_the = klopf_S(Z1, Z2, Gmax, L, f, ereff=ereff)  # Theoretical response
    Sklopf_cas = taperS(Zklopf, L, f, ereff=ereff)         # Actual response of cascade of finite segments

    print(f'Klopfenstein low 3dB cut-off frequency: {klopf_L2f(Z1, Z2, Gmax, L, ereff)*1e-9:.3f} GHz')
    
    fig, axs = plt.subplots(1,2, figsize=(5.5, 5.5/2))
    labels = ['S11', 'S21']
    for inx in [0,1]:
        ax = axs[inx]
        ax.plot(f*1e-9, mag2db(Sklopf_cas[:,inx,0]), '-' , lw=2, label=f'Cascade {N} segments')
        ax.plot(f*1e-9, mag2db(Sklopf_the[:,inx,0]), '--', lw=2, label='Theoretical')
        ax.set_ylabel(labels[inx] + ' (dB)')
        ax.set_xlabel('Frequency (GHz)')
    fig.set_dpi(600)
    fig.tight_layout(pad=0.5)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, bbox_to_anchor=(0.5, 1), 
                   loc='lower center', ncol=2, borderaxespad=0)
    
    # Comparison between taper types...
    # impedance profiles
    Zlin  = (Z2-Z1)*x + Z1    # linear
    Zexp  = Z1*(Z2/Z1)**x     # exponential
    Zatan = (Z2-Z1)*np.arctan(x/(1-x*0.999))*2/np.pi + Z1  # arctan (the 0.999 is to avoid dividing by zero)
    Ztanh = (Z2-Z1)*np.tanh(x/(1-x*0.999)) + Z1            # tanh
    
    fig, ax = plt.subplots(1,1, figsize=(5.5, 5.5/1.5))
    ax.plot(x, Zklopf, '-' , lw=2, label='Klopfenstein')
    ax.plot(x, Zlin  , '-' , lw=2, label='Linear')
    ax.plot(x, Zexp  , '-' , lw=2, label='Exponential')
    ax.plot(x, Zatan , '-' , lw=2, label='Arctan')
    ax.plot(x, Ztanh , '-' , lw=2, label='Tanh')
    ax.set_ylabel('Impedance profile (ohm)')
    ax.set_xlabel('Normalized length')
    fig.set_dpi(600)
    fig.tight_layout(pad=1.2)
    plt.legend()
    
    # frequency response
    Slin  = taperS(Zlin, L, f, ereff=ereff)
    Sexp  = taperS(Zexp, L, f, ereff=ereff)
    Satan = taperS(Zatan, L, f, ereff=ereff)
    Stanh = taperS(Ztanh, L, f, ereff=ereff)
    
    fig, axs = plt.subplots(1,2, figsize=(5.5, 5.5/2))
    labels = ['S11', 'S21']
    for inx in [0,1]:
        ax = axs[inx]
        ax.plot(f*1e-9, mag2db(Sklopf_cas[:,inx,0]), '-' , lw=2, label='Klopfenstein')
        ax.plot(f*1e-9, mag2db(Slin[:,inx,0]) , '-', lw=2, label='Linear')
        ax.plot(f*1e-9, mag2db(Sexp[:,inx,0]) , '-', lw=2, label='Exponential')
        ax.plot(f*1e-9, mag2db(Satan[:,inx,0]), '-', lw=2, label='Arctan')
        ax.plot(f*1e-9, mag2db(Stanh[:,inx,0]), '-', lw=2, label='Tanh')
        ax.set_ylabel(labels[inx] + ' (dB)')
        ax.set_xlabel('Frequency (GHz)')
    fig.set_dpi(600)
    fig.tight_layout(pad=0.5)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, bbox_to_anchor=(0.5, 1), 
                   loc='lower center', ncol=3, borderaxespad=0)
    
    plt.show()
    
    # EOF