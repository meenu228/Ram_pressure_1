# Intera- cluster data
import numpy as np
import cosmology
import numpy as np
from astropy import units as u
from scipy import integrate
from astropy.constants import G
import Galaxy

class Cluster:
    def __init__(self,para):

        #UNITS TO BE FOLLOWED: Mpc Msun and Myr

        cos = cosmology.Cosmology(para)

        #Coma cluster parameters...
        from astropy import units as u
        from astropy.constants import G
        self._G = G.to('Mpc3/(Msun Myr2)')

        self._z = para['z']
        self._omega_m = para['omega_m']
        self._omega_l = para['omega_l']
        self._omega_k = para['omega_k']
        self._h = para['h']
        self._M_200 = para['m_200']*u.Msun                     # IN UNITS OF SOLAR MASS
        self._R_INNER = para['r_inner']*u.Mpc
        #print("self._R_INNER",self._R_INNER)
        self._beta = para['beta']
        self._del1 = para['del1']
        #self._C = para['concentration']
        self._C = ((6.71/((1.0+self._z)**0.44))*((self._M_200*self._h/(2.0e12*u.Msun))**(-0.092))).value

        f_g = 0.8
        O_b = 0.022 / (0.7 ** 2)
        O_m = 0.316

        self._DELTA = para['del1']*f_g * O_b / O_m  # f_uni/(1.0 - f_uni)#del1*f_uni




        #self._DELTA = self._del1*self._f_uni
        self._omega_m_z = (self._omega_m*((1.0+self._z)**3))/(self._omega_m*((1.0+self._z)**3) + self._omega_k*((1.0+self._z)**2) + self._omega_l)
        d = self._omega_m_z - 1.0
        self._Delta_c = 18.0*(np.pi**2) + 82.0*d - 39.0*(d**2)
        self._R_200 = (0.784*((self._M_200*self._h/(1.0e+8*u.Msun))**(1.0/3.0))*((self._omega_m_z*18.0*np.pi**2/(self._omega_m*self._Delta_c))**(1.0/3.0))*(10.0/(self._h*(1.0+self._z)))*u.kpc).to('Mpc')

        self._R_CORE = self._R_200 / (20.0)  # IN UNITS OF MPC

        self._R_VIRIAL = self._R_200 / self._R_CORE
        #print("self._R_VIRIAL  is",self._R_VIRIAL )self._omega_m*self._Delta_c


        self._V200 =  (23.4*((self._M_200*self._h/(1.0e8*u.Msun))**(1.0/3.0))*(((self._omega_m_z*18.0*(np.pi**2))/(self._omega_m*self._Delta_c))**(1.0/6.0))*((self._h*(1.0+self._z)/10)**(1.0/2.0))*(u.km/u.second)).to('Mpc/Myr')

        #print "\n Cluster Parameters: "

        #print "The R200 of cluster: ",self._R_200,self._V200

        H_z = cos.H(1.0 / (1 + cos._z)) * (1 / u.Myr)

        self._delta_c = (200.0 / 3.0) * ((self._C ** 3) / (np.log(self._C + 1.0) - (self._C / (1 + self._C))))


        self._rho_c =  (3.0*(H_z**2)) / (8.0 * np.pi * (G.to('Mpc3/(Msun Myr2)')))
        


    def f1(self,x,R200,c1):
         FUNCT1 = (x**2)/((1+(x)**2)**(3.0*self._beta/2.0))
         return FUNCT1

    def f2(self,x,R200,c1):
        c=c1
        dummy = 0.0001
        #d = self._omega_m_z - 1.0

        #delta_c = (((self._delta_c/(d+1))*self._omega_m)/3.0)*(self._C**3/(np.log(self._C+1.0) - (self._C/(1 + self._C))))

        FUNCT2 = ((x + 0.000001) ** 2) / ((c * (x + 0.000001)) * ((1.0 + c * (x + 0.000001)) ** 2))

        return FUNCT2


    def rho_beta(self,x,IGM,c1):

        R_200= 1.0*IGM._R_200.value
        M_BG = integrate.romberg(self.f2,self._R_INNER.value,1.0,args=(R_200,IGM._C),divmax = 100)
        CONST = integrate.romberg(self.f1,self._R_INNER.value,R_200/self._R_CORE.value,args=(R_200,IGM._C),divmax = 100)
        rho_0 = 0.5*(self._rho_c*self._delta_c)*self._DELTA*(20.0**3)*M_BG/CONST

        rho_beta = (rho_0*((1.0 +(x/self._R_CORE)**2)**(-3.0*self._beta/2.0)))
        #rho_beta  = 277.3362679984032 *(u.Msun *u.kpc**(-2))
        #rho_beta = rho_beta .to (u.Msun *u.Mpc**(-2))
       # print("rho_beta",rho_beta)
        
        

        #print ("rho_beta is",self._rho_c)

        return rho_beta
    
    
  
    def RAM_F(self,x1,v1,IGM):
        Ram_f = self.rho_beta(x1,IGM,IGM._C) * v1 * v1#*np.sin(th) #RPS
        
        return Ram_f
        '''
        M = 10**12 *u.Msun
        print("M",M)
        
        m = 10**11 *u.Msun
        miu = self._G*(M+m)
        print("miu",miu)
        r = 780*u.kpc
        r = r.to(u.Mpc)
        print("r",r)
        R_max = 780*u.kpc
        R_max = R_max.to (u.Mpc)
        print("R_max",R_max)
        e = 0.5
        a = R_max/(1+e)
        print("a",a)
        v1 = (miu*((2/r)-(1/a)))**(0.5) *('u.Mpc*u.Myr**(-1)')
        R = 15 *u.kpc .to (u.Mpc)
        
        self.rho_beta  = 21080.639573797056 *(u.Msun *u.kpc**(-2))
        self.rho_beta = self.rho_beta .to("Msun * Mpc**(-2)")
        #print("v1 is",v1)
        #Ram_f = (self._G *self.rho_beta(x1,IGM,IGM._C)*M)/(np.pi *r_gal) #tidal
        Ram_f = (self.rho_beta(x1,IGM,IGM._C) * v1 * v1)   # for tidal
        print("Ram_f is", Ram_f)
        '''

        '''
         #print("v1 is", v1)
        #th = (np.pi/180)*(np.pi/6)
        M = 10**12 *u.Msun
        #print( "M is", M)
        T_igm = 8.2
        dx = 15 *u.kpc 
        dx = dx .to (u.Mpc)   
        dy = 15 *u.kpc 
        dy = dy .to(u.Mpc)
        r_gal = 15 *u.kpc
        r_gal = r_gal .to(u.Mpc)
        #print("r_gal is ", r_gal)
        lamda_igm = (23*((T_igm/10**8)**2)) *u.kpc 
        lamda_igm = lamda_igm .to (u.Mpc )
        c_igm = 1450 *(u.km *u.s**(-1))
        c_igm = c_igm .to(u.Mpc *u.Myr**(-1))
        #print(dx,dy,r_gal,lamda_igm,c_igm)
        #Ram_f = ((4.3*np.pi*r_gal* self.rho_beta(x1,IGM,IGM._C)* v1*c_igm*lamda_igm)/(np.pi*(dx*dy)))#viscuss
        '''        
