#!/home/storage1/ph14036/miniconda2/bin/python
# main code
import Galaxy
import IGM
import time
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

ti = time.time()
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
#print("r_c200 is",r_c200)
#print("r_c200 is",r_c200)
r_box = 4.0*r_c200

#print("r_box",r_box)
v_c200 = IGM._V200
#print("v_c200",v_c200)

t_total = r_box/v_c200

v=v_c200*para['vel_factor']
#print("v is ",v)
rdg=4*u.kpc.to(u.Mpc)
print("rdg is",rdg)
rds=4*u.kpc.to(u.Mpc)
#print("rds is",rds)
Rout=40*(u.kpc)
Rout = Rout .to (u.Mpc)

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
    #print("acc",acc)
    return A
n = para['n']
n1= para['n1']

outfile = para['outfile']
fp3 = open(outfile,"w")
dm=0.0*u.Msun
dm2 = 0.0*u.Msun
dm3 = 0.0*u.Msun
dm4 = 0.0*u.Msun
dm5 = 0.0*u.Msun
#mapsize1=0
#mapsize3=50
#mapsize4=50
#mapsize2=0

#fp3.close() 
d_i = 0.0*u.Mpc
d_f = 0.0*u.Mpc
u_i = 0.0*u.Mpc/u.Myr
u_f = 0.0*u.Mpc/u.Myr
rc = r_c200
print("rc and v_c200 is",rc,v_c200)
mass_rem1=0.0
mapsize3 = n
mapsize4 = n
dx = 2*Rout/n
dy = 2*Rout/n
#global = my_data

for j in range(1,n1+1):
   # ac=0
    print("j is",j)
    for ii in range(0,n):
        for jj in range(0,n):
            my_data = genfromtxt('sigg1.csv', delimiter=',')
            print(ii,jj)            
            if(gal.sig_g(ii,jj)!=0):
               # th = (np.pi/180)*(np.pi/6)
                if (acc(rc,v,ii,jj) > 0.0):
                   # ac=acc(rc,v,ii,jj)
                    u_f = u_i + acc(abs(rc),v,ii,jj) * dt
                    delta_d = (u_i * dt + 0.50 * acc(abs(rc),v,ii,jj) * (dt ** 2.0))
                    d_f = d_i + delta_d
                    u_i = u_f           
                    d_i = d_f                    
                    if ((gal._z_d - d_i).value > gal._eps):
                        continue
                    else:
                        dm1_1 =  gal.sig_g(ii,jj)#*(u.Msun/u.kpc**2).to('Msun/Mpc**2')
                        dm1_2 =  gal.sig_g (mapsize3-1-jj,ii)
                        dm1_3 =  gal.sig_g (mapsize3-1-ii,mapsize3-jj-1)
                        dm1_4 =  gal.sig_g (jj,mapsize3-1-ii)
                        dm2_1 = dm1_1*((dx*dy))
                        dm2_2 = dm1_2*((dx*dy))
                        dm2_3 = dm1_3*((dx*dy))
                        dm2_4 = dm1_4*((dx*dy))
                        Mdg = 1135724388.0485163*u.Msun

                        dm4 = dm4+(dm2_1/Mdg)*u.Msun
                        dm4 = dm4+(dm2_2/Mdg)*u.Msun
                        dm4 = dm4+(dm2_3/Mdg)*u.Msun
                        dm4 = dm4+(dm2_4/Mdg)*u.Msun
                       # dm5 = dm5+(dm2_1/gal.Mdg+dm2_2/gal.Mdg+dm2_3/gal.Mdg+dm2_4/gal.Mdg)*u.Msun
                       # print(dm2_1+dm2_2+dm2_3+dm2_4)
                       # print("dm4",dm4)
                        
                        im = 0                           
                        my_data[ii,jj]=im
                        my_data[mapsize3-1-jj,ii] = 0
                        my_data[mapsize3-1-ii,mapsize3-jj-1]=0
                        my_data[jj,mapsize3-1-ii]=0
                        #xi = abs((mapsize3-1-ii))
                        #yj = abs((mapsize4-1-jj))
                        #my_data[xi,yj]=im
                        np.savetxt("sigg1.csv",my_data, delimiter=',') #dm.value/gal._m_d.value
                     
                        #print (abs(rc)/IGM._R_200,dm.value/gal._md_tot.value,dm.value,dm.value/gal._m_d.value,file=open('meenu.txt','a'))
                        break
    rc=rc-v*dt
    print(abs(rc)/IGM._R_200,dm4.value,file = open('meenu1.txt','a'))   


    
   # print("rc is",rc,"acc",ac)
    if (j == 20):
        np.savetxt("sig_20.csv",my_data,delimiter = ",")
    if (j==38):
        np.savetxt("sig_38.csv",my_data,delimiter = ",")
    if (j==40):
        np.savetxt("sig_40.csv",my_data,delimiter = ",")    
    if (j == 42):
        np.savetxt("sig_42.csv",my_data,delimiter = ",")
    if (j == 44):
        np.savetxt("sig_44.csv",my_data,delimiter = ",")
    if (j == 46):
        np.savetxt("sig_46.csv",my_data,delimiter = ",")
    if (j == 48):
        np.savetxt("sig_48.csv",my_data,delimiter = ",")            
    if (j == 50):
        np.savetxt("sig_50.csv",my_data,delimiter = ",")
    if (j == 52):
        np.savetxt("sig_52.csv",my_data,delimiter = ",")
    if (j == 55):
        np.savetxt("sig_55.csv",my_data,delimiter = ",")
    if (j == 58):
        np.savetxt("sig_58.csv",my_data,delimiter = ",")
    if (j == 60):
        np.savetxt("sig_60.csv",my_data,delimiter = ",")
    if (j == 62):
        np.savetxt("sig_62.csv",my_data,delimiter = ",")
    if (j == 65):
        np.savetxt("sig_65.csv",my_data,delimiter = ",")
    if (j == 68):
        np.savetxt("sig_68.csv",my_data,delimiter = ",")
    if (j == 70):
        np.savetxt("sig_70.csv",my_data,delimiter = ",")
    if (j == 73):
        np.savetxt("sig_73.csv",my_data,delimiter = ",")
    if (j == 75):
        np.savetxt("sig_75.csv",my_data,delimiter = ",")
    if (j == 78):
        np.savetxt("sig_78.csv",my_data,delimiter = ",")
    if (j == 80):
        np.savetxt("sig_80.csv",my_data,delimiter = ",")
    if (j == 88):
        np.savetxt("sig_88.csv",my_data,delimiter = ",")
    if (j ==99):
        np.savetxt("sig_99.csv",my_data,delimiter = ",")
    if (j == 100):
        np.savetxt("sig_100.csv",my_data,delimiter = ",")
    if (j == 112):
        np.savetxt("sig_112.csv",my_data,delimiter = ",")
    if (j == 125):
        np.savetxt("sig_125.csv",my_data,delimiter = ",")
    if (j == 140):
        np.savetxt("sig_140.csv",my_data,delimiter = ",")    
    
    
elapsed = time.time()-ti
print(elapsed)       
#fp3.close()
 
