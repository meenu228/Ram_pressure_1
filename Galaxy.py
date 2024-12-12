#galaxy data
import numpy as np
from scipy import integrate
from astropy import units as u
import cosmology
import pandas as pd 

class Galaxy:
    def __init__(self, para):
        from astropy.constants import G

        cos = cosmology.Cosmology(para)

        self._G = G.to('Mpc3/(Msun Myr2)')
        #print(""self._G)

        self._lmbda = para['lmbda']
        self._z = para['z']
        self._alpha = para['alpha']
        self._f_g = para['f_g']
        self._md_tot = para['md_tot'] * u.Msun
        print(self._md_tot)
        self._z_d = para['z_d'] * u.Mpc
        self._m_bh = para['m_bh'] * u.Msun
        self._m_bul = para['m_bul'] * u.Msun
        self._m_tot = self._m_bh + self._m_bul
        self._epsilon = para['epsilon']
        self._eps = para['eps']
        self._f_uni = para['f_uni']
        self._z = para['z']
        self._omega_m = para['omega_m']
        self._omega_l = para['omega_l']
        self._omega_k = para['omega_k']
        self._h = para['h']
        #print("self._eps",self._eps)
        self._omega_m_z = (self._omega_m * ((1.0 + self._z) ** 3)) / (
        self._omega_m * ((1.0 + self._z) ** 3) + self._omega_k * ((1.0 + self._z) ** 2) + self._omega_l)
        d = self._omega_m_z - 1.0
        self._delta_c = 18.0 * np.pi ** 2 + 82.0 * d - 39.0 * d ** 2

        # G_gal = G.to('pc3/Msun yr2')

        self._m_d = self._f_g * self._f_uni * self._md_tot
        print("self._m_d is",self._m_d)

        #print 'Total mass is galaxy is', self._md_tot, "total mass in disk = ", self._alpha * self._f_uni * self._md_tot, "mass is gas= ", self._f_g * self._m_d, "Mass in stars: ", (1.0 - self._f_g) * self._m_d
        self._r_dg = ((self._lmbda/np.sqrt(2.0))*0.784*((self._md_tot*self._h/(1.0e+8*u.Msun))**(1.0/3.0))*((self._omega_m_z*18.0*np.pi**2/(self._omega_m*self._delta_c))**(1.0/3.0))*(10.0/(self._h*(1.0+self._z)))*u.kpc).to('Mpc')

        H = (cos.H(1.0 / (1 + cos._z))) * (1 / u.Myr)

        #self._r_dg = (self._lmbda / np.sqrt(2.0)) * ((self._md_tot * 10.0 * self._G * H) ** (1.0 / 3.0)) / (10.0 * H)

        self._r_out = 10.0 * self._r_dg

        self._r_ds = self._r_dg / 2.0  # 3.0*self._r_out/20.0
        self.Mdg = self._f_g *self._m_d
        #print("self.Mdg is",self.Mdg)
        self.Mds = (1.0-self._f_g)*self._m_d
        #print("self.Mds is",self.Mds)

    #def funct_sig_g(self, y1):
     #   f1 = y1 * np.exp(-y1)
      #  return f1

    #def funct_sig_s(self, y2):
     #   f2 = y2 * np.exp(-y2)
      #  return f2

    def sig_g(self, ii,jj):
        Gas_data = np.loadtxt('sigg1.csv',delimiter=',')
        Gas_data = (Gas_data*(u.Msun/u.kpc**2)).to('Msun/Mpc**2')
        #print(Gas_data)
        #print(Gas_data[ii,jj])
        return Gas_data[ii,jj]

    def sig_s(self, ii,jj):
        Star_data = np.loadtxt('sigs1.csv',delimiter=',')
        #Sig_s = Star_data*u.Msun/((0.4*(u.kpc)).to('Mpc'))**2
        Sig_s = (Star_data*(u.Msun/u.kpc**2)).to('Msun/Mpc**2')
        #Sig_s = Sig_s*u.Msun/((u.kpc).to('Mpc'))**2 #conversion factor from pixels
        
        return Sig_s[ii,jj]
        #return star_data[ii,jj]
    
    def Magnetic_F(self,ii,jj):
        Force = np.loadtxt("mag_s.csv",delimiter = ',')
        Force = Force*u.Msun /((u.Mpc)*(u.Myr**2))
        #print("Lorent Force shape is",Lorent_Force.shape)
        return Force[ii,jj]
       
    def RES_F(self, ii,jj):
        #print "M here",self.sig_g(ii,jj),self.sig_s(ii,jj)
        f_s = (2.0 * np.pi * (self._G) * self.sig_g(ii,jj) * self.sig_s(ii,jj)) + self.Magnetic_F(ii,jj)
        #print(f_s)
        #f_b = (self._G) * self._m_tot * self.sig_g(ii,jj) / (x2 ** 2.0)
        Res_f = (f_s + self._epsilon)
        #print(Res_f)
        return Res_f

