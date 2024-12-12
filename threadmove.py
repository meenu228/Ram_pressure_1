#!/home/storage1/ph14036/miniconda2/bin/python
import Galaxy
import IGM
from configparser import SafeConfigParser
from astropy import units as u
import numpy as np
from numpy import genfromtxt
from scipy import integrate
from numba import jit
import sys
from scipy.integrate import quad
from scipy.integrate import dblquad
from math import *
import matplotlib.pyplot as plt
import math
import csv
import pandas as pd
import threading

parser = SafeConfigParser()
parser.read(sys.argv[1])

para = {}

from astropy.constants import G

G = G.to('Mpc3/(Msun Myr2)')

#print('\nCosmology:')
for name in parser.options('cosmology'):
    string_value = parser.get('cosmology', name)
    value = parser.getfloat('cosmology', name)
    para[name] = value
   # print('  %-12s : %-7r -> %0.2f' % (name, string_value, value))

#print('\nCluster')
for name in parser.options('cluster'):
    string_value = parser.get('cluster', name)
    value = parser.getfloat('cluster', name)
    para[name] = value
  #  print('  %-12s : %-7r -> %0.2f' % (name, string_value, value))

##print('\ngalaxy')
for name in parser.options('galaxy'):
    string_value = parser.get('galaxy', name)
    value = parser.getfloat('galaxy', name)
    para[name] = value
  #  print('  %-12s : %-7r -> %0.2f' % (name, string_value, value))

#print('\nmove')
for name in parser.options('move'):
    string_value = parser.get('move', name)
    if name=='n1':
        value = parser.getint('move', name)
    elif name=='n':
        value = parser.getint('move', name)
    elif name =='z':
        value = parser.getfloat('move', name)
    else:
        value = parser.get('move', name)

    para[name] = value
    #print '  %-12s : %-7r -> %0.2f' % (name, string_value, value)

#print("n1,n = ",para['n1'],para['n'])

gal = Galaxy.Galaxy(para)
IGM = IGM.Cluster(para)

r_c200 = IGM._R_200
r_box = 4.0*r_c200
v_c200 = IGM._V200
t_total = r_box/v_c200

v=v_c200*para['vel_factor']
rdg=2.5*u.kpc.to(u.Mpc)
rds=2.5*u.kpc.to(u.Mpc)
Rout=15
mass=gal._m_d/u.Msun
#IGM._R_VIRIAL = IGM._R_VIRIAL*u.Mpc
#print("IGM._R_VIRIAL",IGM._R_VIRIAL)
#print("total time taken to cross cluster: ",(r_box/v_c200).to('Myr'))
#print("es is",es)
dt = t_total/para['n1']
dr = gal._r_out/para['n']
#####################################################
########################################################

def acc(x1,v1,ii,jj):
    #print IGM.RAM_F(x1,v1,IGM),gal.RES_F(ii,jj)
    A = (IGM.RAM_F(x1,v1,IGM) - gal.RES_F(ii,jj))/(gal.sig_g(ii,jj))
    #print IGM.RAM_F(x1,v1),gal.RES_F(x2),gal.sig_g(x2)
    return A
n = para['n']
n1= para['n1']
outfile = para['outfile']
fp3 = open(outfile,"w")
dm=0.0*u.Msun
map_size= 50
fp3.close() 
d_i = 0.0*u.Mpc
d_f = 0.0*u.Mpc
u_i = 0.0*u.Mpc/u.Myr
u_f = 0.0*u.Mpc/u.Myr
rc = r_c200
#mass_rem1=0.0
dx = 2*Rout/100
dy = 2*Rout/100

def thread(ii,jj):
    global my_data
    global u_i
    global d_f
    global d_i
    global u_f
    
    #th = (np.pi/180)*(np.pi/6)
    if(gal.sig_g(ii,jj)!=0):
        #th = (np.pi/180)*(np.pi/6)
        if (acc(rc,v,ii,jj) > 0.0):
            u_f = u_i + acc(abs(rc),v,ii,jj) * dt
            delta_d = (u_i * dt + 0.50 * acc(abs(rc),v,ii,jj) * (dt ** 2.0))
            d_f = d_i + delta_d
            u_i = u_f           
            d_i = d_f
            if not((gal._z_d - d_i).value > gal._eps):
                
            
                my_data = genfromtxt('sigg1.csv', delimiter=',')
                dm =  gal.sig_g(ii,jj)*((0.075*(u.kpc)).to('Mpc'))**2
                dm = dm/(dx*dy)
                #print("dm is",dm)
                #sig = gal.sig_g(ii,jj)                      
                im=my_data[ii,jj]
                #print("im is",im)
                if(im-dm.value<0):
                    my_data[ii,jj]=0
                else:
                    my_data[ii,jj]=im-dm.value
               # my_data[ii,jj]=im-dm.value
                np.savetxt("sigg1.csv",my_data, delimiter=',')
                #print(my_data[ii,jj])
                
    
my_data = genfromtxt('sigg1.csv', delimiter=',')  
for j in range(1,100):
    print("j is",j)
    for ii in range(0,100):
        for jj in range(0,100):
            t=threading.Thread(target=thread, args=(ii,jj))
            t.start()
            t.join()
            print(ii,jj)
            
    if (j == 35):
        np.savetxt("sigg_35.csv",my_data,delimiter = ",")
    if (j == 65):
        np.savetxt("sig_65.csv",my_data,delimiter = ",")
    if (j == 99):
        np.savetxt("sig_99.csv",my_data,delimiter = ",")
   # if (j == 40):
    #    np.savetxt("sig_50.csv",my_data,delimiter = ",")     
    rc=rc-v*dt
                
           

                    #sig = sig-dm
                    # Rstrip= np.sqrt((-rout+(ii*dr1))**2+(-rout+(jj*dr1))**2)
                    # f8 = lambda x: (x*(gal.sig_g(ii,jj)*u.Mpc**2/u.Msun))		 
                    # res_m=integrate.quad(f8,Rstrip,rout)
                    # res_m2=res_m[0]
                    #mass_rem=2*np.pi*res_m2/(0.1*mass)
                    #mass_rem1=mass_rem1+mass_rem
                    # print(("{}\t{}".format(abs(rc)/IGM._R_200,Rstrip/rdg)),file=open('M1.txt','a'))
