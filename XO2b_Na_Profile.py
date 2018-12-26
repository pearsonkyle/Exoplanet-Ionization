import pickle 
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import brentq as findzero

from Ionization import ion, NaI

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
scaleheight = lambda T,u,g : (kb*T)/(u*mp*g) # u - mean mass [amu]
gravity = lambda M,r : G*M/r**2
altitude = lambda H, Pz, Po :  -H * np.log(Pz/Po)
density = lambda P, T : (P*1e6) / (kb*T) # pressure in bars -> barye
ev2nm = lambda e: 1239.84193/e
nm2ev = lambda n: 1239.84193/n

class interpolate():
    def __init__(self,data):
        self.data = np.array(data)
        self.fn = interp1d(self.data[:,0],self.data[:,1],assume_sorted=False,kind='linear',fill_value='extrapolate')
    def __call__(self,P):
        return self.fn(P)
    @property
    def x(self):
        return self.data[:,0]
    @property
    def y(self):
        return self.data[:,1]

class planet():
    '''
    grids
        pressure
        altitude
        gravity
        temperature
        scaleheight
        density
        optical_depth 
        photoionization 
        thermal_ionization
        balance
        recombination
    '''

    opacity = []

    def __init__(self,pars=None, TP=None):
        self.pars = pars
        self.TP = TP 

    def initialize_1d(self):

        # compute uniform spacing in log10(pressure) scale 
        Pmin = self.TP.data[:,0].min()
        Pmax = self.TP.data[:,0].max()
        Hest = kb*self.TP(Pmax) / (self.pars['mu']*mp*gravity(self.pars['mass'],self.pars['radius']) )
        zmax = 1.5 * -Hest * np.log(Pmin/Pmax)

        self.pressure = np.logspace(np.log10(Pmax),np.log10(Pmin),self.pars.nlayers)
        self.altitude = np.zeros(self.pars.nlayers)

        # find altitudes for each pressure level 
        for i in range(self.pressure.shape[0]):
            Pz = self.pressure[i]
            Po = Pmax # pressure at z=0, 10 bar pressure usually
            Tz = self.TP(Pz)
            
            # find alitude with varying P, TP Profile & Gravity 
            # find roots to: T(Pz)*[ln(Pz) - ln(Po)] = -z mu mp * GM / (kb*(R+z)^2) 

            def f2zero(z):
                side1 = Tz*(np.log(Pz) - np.log(Po)) 
                side2 = -z * self.pars.mu * mp * G * self.pars.mass / (kb*(self.pars.radius+z)**2)
                return side1-side2
            
            root = findzero(f2zero, 0, zmax)
            self.altitude[i] = root 
            
        # check altitude with constant scale height and gravity approx below 
        #const = -Hest * np.log(self.pressure /Pmax)

        self.gravity = gravity(self.pars.mass,self.altitude+self.pars.radius)
        self.temperature = self.TP(self.pressure)
        self.scaleheight = scaleheight(self.temperature, self.pars.mu, self.gravity)
        self.density = density(self.pressure,self.temperature)

class AttributeDict(dict):
    def __getattr__(self, name):
        return self[name]
    def __setattr__(self,name,val):
        self[name] = val

if __name__ == "__main__":

    # parameters for planet 
    pars = AttributeDict({
        'radius':0.97*rjup,
        'distance':0.0386, # AU
        'mass':0.567*mjup, 
        'mu':2.3,
        'nlayers':100,
    })

    # create pressure grid and get temperature
    # Pressure(bar), Temperature (K)
    PT = interpolate([ 
        [10,1800],
        [0.005, 1200],
        [1e-10, 1200],
    ])

    # load data from star into grid 
    stellar = np.loadtxt("full_solar.dat")
    binsize = stellar[:,1]+1 - stellar[:,0]
    avgwave = (stellar[:,1]+1 + stellar[:,0])*0.5
    flux = interpolate( np.array([avgwave, (1./pars.distance)**2 * stellar[:,2]/binsize]).T ) # N/cm^2/s/nm

    # ion limit = 241.2 nm
    nm = np.linspace(1,241,121)
    ev = nm2ev(nm)

    # initialize planet grid
    XO2b = planet(pars,PT)
    XO2b.initialize_1d()

    mr = 4.08e-7 # starting mixing ratio
    ni = np.zeros(XO2b.pars.nlayers)

    # create new grids
    XO2b.optical_depth = np.zeros(XO2b.pars.nlayers)
    XO2b.photoionization = np.zeros(XO2b.pars.nlayers)
    XO2b.thermal_ionization = np.zeros(XO2b.pars.nlayers)
    XO2b.balance = np.zeros(XO2b.pars.nlayers)
    XO2b.recombination = np.zeros(XO2b.pars.nlayers)

    # production rate of ions
    def photorate(opt_ni, *args):
        
        # set electron abundance to value from optimization function
        opt_layer, = args
        ni[opt_layer] = opt_ni

        # compute neutral sodium abundance
        nna = mr - ni

        # compute photoionization but first compute optical depth across ionizing wavelengths
        tau = np.zeros( nm.shape[0] ) 
        cross_sections = ion.ionization_crosssection(ev, NaI) # cm^2
        for i in range(nm.shape[0]):
            tau[i] = np.trapz( (nna*XO2b.density*cross_sections[i])[opt_layer:], XO2b.altitude[opt_layer:] ) 
            
        phot_ion, phot_rate = ion.ionization_photo(nm, flux, tau, (nna*XO2b.density)[opt_layer], cross_sections) # N/s/cm3
        
        # rates for each layer in the atmosphere
        recom_rad =  ion.recombination_radiative(XO2b.temperature, NaI)
        recom_3body = ion.recombination_3body(XO2b.temperature, XO2b.density[opt_layer], NaI)
        thermal_ion = ion.ionization_thermal(XO2b.temperature, NaI)

        # total ion production rate, equilibrium = 0
        recom_rate = recom_rad * (ni*XO2b.density)**2 
        thermal_rate = thermal_ion * (ni*XO2b.density)*(nna*XO2b.density)
        body3_rate = recom_3body * (ni*XO2b.density)**2
        balance = phot_ion - recom_rate + thermal_rate - body3_rate

        # save values to planetary atmo grid
        XO2b.balance[opt_layer] = balance[opt_layer]
        XO2b.photoionization[opt_layer] = phot_rate
        XO2b.recombination[opt_layer] = recom_rate[opt_layer]
        XO2b.thermal_ionization[opt_layer] = thermal_rate[opt_layer]

        # return balance at specific layer 
        return balance[opt_layer]

    # optimize the mixing ratio in each atmospheric layer
    # start at the bottom of the atmosphere and work our way up
    for layer in range( XO2b.pars.nlayers ):

        ni[layer] = findzero( photorate, 0, mr, args=(layer,) )
        print("{:.2f} log10(bar) {:.2e} [Na+]".format(np.log10(XO2b.pressure[layer]),
                                                            ni[layer]) )

        # update XO2b grids with final value
        photorate(ni[layer], layer)

    # make some plots
    LOW, UPPER = 10**1, 10**-7
    f,ax = plt.subplots(1,3, figsize=(9,4))
    plt.subplots_adjust(left=0.10,right=0.91,bottom=0.19,top=0.84, wspace=0.12)

    ax[0].semilogy(XO2b.temperature,XO2b.pressure,'k-') 
    ax[0].set_ylim([LOW,UPPER])
    ax[0].set_ylabel(r"Pressure [bar]")
    ax[0].set_xlabel("Temperature [K]")
    nna = mr - ni

    ax[1].loglog(ni,XO2b.pressure,'k--',label=r'Na$^{+}$') 
    ax[1].loglog(nna,XO2b.pressure,'k-',label='Na') 
    ax[1].set_ylim([LOW,UPPER])
    ax[1].set_xlim([1e-11,5e-7])
    ax[1].yaxis.set_ticklabels([])
    ax[1].set_xlabel("Mixing Ratio")
    ax[1].legend(loc='best')

    ax[2].loglog(XO2b.photoionization,XO2b.pressure,'k-',label='Photo') 
    ax[2].loglog(XO2b.thermal_ionization,XO2b.pressure,'k--',label='Thermal') 
    ax[2].set_ylim([LOW,UPPER])
    ax[2].set_xlim([10e-10,10e-2])
    ax[2].legend(loc='best')
    ax[2].set_xlabel(r"Ionization Rate [s$^{-1}$]")
    ax[2].set_ylabel(r"Pressure [bar]")
    ax[2].yaxis.tick_right()
    ax[2].yaxis.set_label_position("right")
    plt.show()

    '''
    # atmospheric optical depth
    waves = np.linspace(50,200,3)
    ev = nm2ev( waves )
    cross_sections = ion.ion_cross_section(ev, ion.NaI) # m^2

    for j in range(len(cross_sections)):

        tau = np.zeros(XO2b.pars.nlayers)
        for i in range(XO2b.pars.nlayers):
                tau[i] = trapz( (nna*XO2b.density*cross_sections[j])[i:], XO2b.altitude[i:] ) 
        
        # TODO integrate upwards 
        ax[3].semilogy(np.exp(-tau),XO2b.pressure,label='{:.1f} nm'.format(waves[j]) ) 
        ax[3].set_ylim([10*1,10**-10])
        ax[3].set_xlabel("Atmospheric Opacity")
        ax[3].legend(loc='best')
    
    ax[3].semilogy(np.exp(-1)*np.ones(XO2b.pars.nlayers),XO2b.pressure,'k--' )
    '''

