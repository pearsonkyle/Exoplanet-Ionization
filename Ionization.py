from functools import wraps 
import numpy as np

kb = 1.38e-16  # g cm2 s-2 K-1
G = 6.674e-8 # cm3 g-1 s-2
h = 4.135667516e-15 # eV s
c = 2.9979e10 # cm s-1
# barye = g cm-1 s-2
rsun = 6.957e10 # cm 
rjup = 6.9911e9 # cm 
mjup = 1.898e30 # g
mp = 1.67e-24 # g

# basic functions
density = lambda P, T : (P*1e6) / (kb*T) # pressure in bars -> barye
ev2nm = lambda e: 1239.84193/e
nm2ev = lambda n: 1239.84193/n

# If you would like additional atoms please follow the same format as the 
# dictionaries below. The keys are case and probably order sensitive, so keep everything
# the same format. Coefficients can be found in the following references:
# http://adsabs.harvard.edu/abs/1996ApJS..103..467V
# http://adsabs.harvard.edu/abs/1996ApJ...465..487V
# http://adsabs.harvard.edu/abs/1997ADNDT..65....1V
# See Readme.md for more details 

NaI = {
    'ionization_crosssection':{'E0':6.139,'sigma0':1.601,'ya':6.148e3,'P':3.839,'yw':0,'y0':0,'y1':0}, # cm^2
    'recombination_radiative':{'a':5.641e-12,'b':0.1749,'T0':3.077e2,'T1':2.617e6}, # cm^3 s^-1
    'recombination_3body':{'k1':1e-7,'k0':lambda T:3.43e-14*T**-3.77 },   # cm^3 s^-1
    'ionization_thermal':{'dE':5.1,'P':1,'A':0.101e-6,'X':0.275,'k':0.23} # cm^3/s
}

HI = {
    'ionization_crosssection':{'E0':4.298e-1,'sigma0':5.475e4,'ya':3.288e1,'P':2.963,'yw':0,'y0':0,'y1':0}
    #'ionization_crosssection':{'E0':6.139,'sigma0':1.601,'ya':6.148e3,'P':3.839,'yw':0,'y0':0,'y1':0},
}


class ion():
    '''
    A fancy class that is more like container for functions related to ionization
    These functions can be used in the traditional sense or passed a dictionary like the ones above

    # able to call a function both ways
        ewaves = nm2ev( np.linspace(5,121,40) )
        cross = ion.ionization_crosssection( ewaves, HI)
        cross = ion.ionization_crosssection( ewaves, 4.298e-1,5.475e4,3.288e1,2.963,0,0,0)
    '''

    # call function corresponding to key name with children keys as kw arguments
    def _get_coefficients(foo):
        @wraps(foo)
        def magic( *args ):
            if isinstance(args[-1],dict):
                pars = args[-1].get(foo.__name__)
                return foo(*args[:-1],**pars)
            else:
                return foo(*args)
        return magic

    @staticmethod
    @_get_coefficients
    def ionization_crosssection(E,  E0,sigma0,ya,P,yw,y0,y1):
        # Verner 1996b
        # returns cross section in cm^2
        x = E/E0 - y0
        y = np.sqrt(x**2 + y1**2)
        F = ((x-1)**2 + yw**2) * (y**(0.5*P-5.5)) * (1+np.sqrt(y/ya))**-P
        return sigma0*F*10**-18 

    @staticmethod
    @_get_coefficients
    def recombination_radiative(T,  T0,T1,a,b):
        # Verner 1996a
        # returns recomb rate cm^3/s
        tito = np.sqrt(T/T0)
        ar = a*(tito*(1+tito)**(1-b) * (1+np.sqrt(T/T1))**(1+b)  )**-1 
        return ar 

    @staticmethod
    @_get_coefficients
    def recombination_3body(T,ni,  k1,k0):
        #  where ni is ion density in cm-3
        # 1 = 1.d-7, k0 = (3.43d-14)*(temp^-3.77) coefficients are for cgs units
        return ( k1*k0(T)*ni/(k0(T)*ni + k1) )

    @staticmethod
    @_get_coefficients
    def ionization_thermal(T,  dE,P,A,X,k):
        # Voronov 1997
        # terminal ionization rate by electron impact, cm^3/s
        
        # convert temperature [K] to eV
        TE = T*8.617343e-5
        U = dE/TE
        ov = A * (1 + P*U**0.5) * U**k * np.exp(-U) / (X+U)
        return ov 

    @staticmethod
    def ionization_photo(wave,flux,tau,dens,cross):    
        # integrates the photoionization rate over wavelengths 
        # wave - nm
        # flux - phot/cm2/s
        # density - cm-3
        # cross section - cm2
        photo = np.zeros( wave.shape[0])
        phot_rate = np.zeros( wave.shape[0])
    
        for i in range(len(wave)):
            # photons/cm2/s from star
            f = flux( wave[i] ) 

            # ionization rate [N/s/cm^3] per wavelength
            photo[i] = f * np.exp(-tau[i]) * dens * cross[i]

            # ionization rate [N/s] per wavelength
            phot_rate[i] = f * np.exp(-tau[i]) * cross[i]

        # integrate rate over all wavelengths
        return np.trapz( photo, wave ), np.trapz( phot_rate, wave )



if __name__ == "__main__":

    # able to call a function both ways
    #ewaves = nm2ev( np.linspace(5,121,40) )
    #cross = ion.ionization_crosssection( ewaves, HI)
    #cross = ion.ionization_crosssection( ewaves, 4.298e-1,5.475e4,3.288e1,2.963,0,0,0)

    # define our atmospheric parameters 
    T, P = 1500, 10  # [K], [bar]
    ni = 1e-10       # Na ion mixing ratio (same as electron mr)
    nna = 3e-6 - ni  # neutral Na mixing ratio 
 
    # function will automatically get coefficients from ion.NaI dictionary 
    thermal_rate = ion.ionization_thermal(T, NaI) # cm^3 / s 

    # rate coefficient * neutral Na density * electron density 
    thermal_ionization = thermal_rate * (nna*density(P,T)) * (ni*density(P,T))

    # compute radiative recombination rate
    recom_rate =  ion.recombination_radiative(T, NaI)

    # assume ion density = electron density 
    radiative_recombination = recom_rate * (ni*density(P,T))**2

    print("thermal rate coefficient <ov> [cm^3 s^-1] {:.2e}".format(thermal_rate))
    print("thermal ionization production [s^-1 cm^-3] {:.2e}\n".format(thermal_ionization) )

    print("rad recombination rate coeff  [cm^3 s^-1] {:.2e}".format(recom_rate))
    print("radiative recombination prod  [s^-1 cm^-3] {:.2e}".format(radiative_recombination) )