
# coding: utf-8

# # Poincare section

# Importing libraries, declaring constants

# In[ ]:

import numpy as np
import scipy as sp
import math
from math import atan, sin, cos, acosh, exp, log, cosh, tanh, pi
import os
from os import mkdir, makedirs

delta = 1.0
eps = 1e-7 #for energy conservation in Runge-Kutta method


# Defining field components

# In[ ]:

def Field(x, y, z):
         
    field_components = dict()
    
    field_components['bx'] = bx0 * tanh(z / delta)  
    field_components['bz'] = bn
    
    if fieldtype == 'constant' or fieldtype == 'shearless':
        field_components['by'] = by0
        
    elif fieldtype == 'bell':
        if abs(z) < delta:
            field_components['by'] = by0 * cos( (z / delta) * (pi / 2.0))
        else:
            field_components['by'] = 0.0
    
    elif fieldtype == 'antisymm':
        if abs(z) < delta:
            field_components['by'] = - by0 * sin( z * pi / delta)
        else:
            field_components['by'] = 0.0

    return field_components


# Defining derivative functions for Runge-Kutta method

# In[ ]:

def f(b, a):
        
    q = np.empty((6,))
    q[0] = a[3]
    q[1] = a[4]
    q[2] = a[5]
    q[3] = delta * (b['bz']*a[4] - b['by']*a[5])
    q[4] = delta * (b['bx']*a[5] - b['bz']*a[3])
    q[5] = delta * (b['by']*a[3] - b['bx']*a[4]) 
    
    return q


# Defining Runge-Kutta method

# In[ ]:

def RK4_Poinc(f, a0):
    Cy = a0[4]  + delta*bn*a0[0] - delta**2 * bx0 * log( cosh( a0[2] / delta ))
    a0[0] = a0[0] - Cy / (bn*delta) #shift of coordinates    
    
    
    RK = a0.reshape((1, 6)) #trajectory array
    times = [0.0] #array with times
    Poinc = np.zeros((pcountmax, 2)) #array with found poincare points for the current particle
    ZeroIndices = [] #indices of poincare points in trajectory array
    
    en0 = (a0[3]**2 + a0[4]**2 + a0[5]**2)**0.5 #initial energy equals 1.0
    
    t = 0.0 #time
    i = 0 #current trajectory point
    idiv = 0 #number of times division dt/2 was used        
    pcount = 0 #number of found points suitable for Poincare section
    toofar = False #whether the particle is too far from the layer
    
    OldRK = a0    
       
    while t < tfinal:
        
        ##==============Runge-Kutta step===============
        
        dt = dt0
        for s in range(100):  
            
            b = Field(OldRK[0], OldRK[1], OldRK[2])
            k1 = dt * f(b, OldRK)
            
            b = Field(OldRK[0] + 0.5*k1[0], OldRK[1] + 0.5*k1[1], OldRK[2] + 0.5*k1[2])
            k2 = dt * f(b, OldRK + 0.5*k1)
            
            b = Field(OldRK[0] + 0.5*k2[0], OldRK[1] + 0.5*k2[1], OldRK[2] + 0.5*k2[2])
            k3 = dt * f(b, OldRK + 0.5*k2) 
            
            b = Field(OldRK[0] + k3[0], OldRK[1] + k3[1], OldRK[2] + k3[2])
            k4 = dt * f(b, OldRK + k3) 
            
            NewRK = OldRK + ( k1 + 2.0 * ( k2 + k3 ) + k4 ) / 6.0
            en = (NewRK[3]**2 + NewRK[4]**2 + NewRK[5]**2)**0.5
            
            if abs(en - en0) < eps: #energy conservation control
                break  
            dt = dt / 2.0
            idiv +=1
            
        RK = np.vstack((RK, NewRK)) #adding new line to trajectory array
        t +=dt
        times.append(t)
            
        ##=================== Poincare==================
        
        if OldRK[2] * NewRK[2] < 0: #z coordinate sign change condition
            
            ZeroIndices.append(i)
            
            if pcount < pcountmax:        
                                    
                #linear approximation of bn*delta*(x)
                Poinc[pcount, 0] = bn * delta * (OldRK[0] + OldRK[2] * (NewRK[0] - OldRK[0]) /                                                  (OldRK[2] - NewRK[2]))
                
                #linear approximation of Vx
                Poinc[pcount, 1] = OldRK[3] + OldRK[2] * (NewRK[3] - OldRK[3]) / (OldRK[2] - NewRK[2])
                                    
                pcount +=1
        
        ##============cycle break conditions============
        
        if pcount == pcountmax: #enough poincare points
            break
      
        if abs(NewRK[2]) > toofarcondition: #particle too far from the layer
            toofar = True
            break
              
        ##==============================================
        
        OldRK = NewRK
        i +=1
    
    length = len(times)
    times = np.asarray(times).reshape((length, 1))
    RKt = np.hstack((RK, times)) #array with Runge-Kutta results and times
    
    print 'p =', p
    print 'Initial approximation (x0, y0, z0, Vx0, Vy0, Vz0) = (', round(a0[0], 5), ',', round(a0[1], 5), ',', round(a0[2], 5),           ',', round(a0[3], 5), ',', round(a0[4], 5), ',', round(a0[5], 5), ')'
    print 'Particle too far from the layer =', toofar
    print 'Initial energy =', en0, ', final energy =', en
    print 'Time =', t, ', Total calculated points =', i+1
    print 'Total Poincare points =', len(ZeroIndices), '/', pcountmax
    #print 'Number of times division dt/2 was used =', idiv
    print '\n'
   
    return RKt, Poinc, ZeroIndices 


# Setting

# In[ ]:

root = r'Res2016-07-26'

# In[ ]:

dt0 = 0.01
tfinal = 3000.0
toofarcondition = 12.0 #condition of whether the particle is too far from the layer
pcountmax = 75 # maximum Poincare points
maxpos = 60 #number of different positions +1, must be even
maxspeed = maxpos #number of different initial Vx speeds +1, must be equal to maxpos
h = 2.0 / maxspeed
pmax = int(maxpos * maxspeed * pi / 4)


# ## One particle

# In[ ]:
'''
fieldtype = 'constant'
#fieldtype = 'bell'
#fieldtype = 'antisymm'
#fieldtype = 'shearless'
bx0 = 1    
by0 = 0.5
bn = 0.10

folder = r'%s\%s\by=%.2fbx\bn=%.2fbx\Poincare' %(root, fieldtype, by0, bn)
os.makedirs(folder)

norm = (bx0**2 + by0**2 + bn**2)**0.5
bx0 = bx0 / norm
by0 = by0 / norm
bn = bn / norm

logfile1 = open(r'%s\parameters.txt' %(folder), 'w')
logfile1.write('Field type = %s\nbx0 = %g\nby0 = %g\nbn = %g\n\ndelta = %g\ndt0 = %g\ntfinal = %g\npcountmax = %g\n\n'               %(fieldtype, bx0, by0, bn, delta, dt0, tfinal, pcountmax))
logfile1.write('maxpos = %g\nmaxspeed = %g\n'               %(maxpos, maxspeed))
logfile1.write('number of particles, approximately = %g\n\n' %(pmax)) 
logfile1.close() 

pos = 0
speed = 0
p = 0

a00 = ( -1.0 + 2.0 * (pos * 1.0 / maxpos) )/ (bn * delta) #chosen x0*delta*bn from -1 exlusive to 1 exclusive
a01 = 0.0
a02 = 0.0   
      
a03 = (-1.0 + 2.0 * (speed * 1.0 / maxspeed)) #distribution of Vx                                                                        
a04 = - bn * delta * a00
a05 = (1.0 - a03**2 - a04**2)**0.5 #so that the initial energy equals 1.0, nesessarily > 0
                        
a0 = np.array([a00, a01, a02, a03, a04, a05])

print 'pos =', pos, ', speed =', speed
RKt, Poinc, ZeroIndices = RK4_Poinc(f, a0)            
 
    
np.savetxt(r'%s\Trajectory_%d.txt' %(folder, p),                    RKt, ['%.4g','%.4g','%.4g','%.4g','%.4g','%.4g','%g',], delimiter=',') 
        
indicesfile = open(r'%s\ZeroIndices_%d.txt' %(folder, p), 'w')
for i in ZeroIndices: 
    indicesfile.write('%s\n' %i)    
indicesfile.close()     
        
np.savetxt(r'%s\Poinc_%d.txt' %(folder, p), Poinc, '%.3g', delimiter=',')

'''
# ## Multiple particles (for Poincare section)

# In[ ]:

for fieldtype in ['shearless', 'constant', 'bell', 'antisymm']:
    if fieldtype == 'shearless':
        by0s = [0.0]
    else:
        by0s = [0.25, 0.5, 0.75, 1]
    for by0 in by0s:
        theby0 = by0
        for bn in [0.1, 0.2]:
            
            bx0 = 1   
            by0 = theby0 

            folder = r'%s\%s\by=%.2fbx\bn=%.2fbx\Poincare' %(root, fieldtype, by0, bn)
            os.makedirs(folder)

            norm = (bx0**2 + by0**2 + bn**2)**0.5
            bx0 = bx0 / norm
            by0 = by0 / norm
            bn = bn / norm

            logfile1 = open(r'%s\parameters.txt' %(folder), 'w')
            logfile1.write('Field type = %s\nbx0 = %g\nby0 = %g\nbn = %g\n\ndelta = %g\ndt0 = %g\ntfinal = %g\npcountmax = %g\n\n'                           %(fieldtype, bx0, by0, bn, delta, dt0, tfinal, pcountmax))
            logfile1.write('maxpos = %g\nmaxspeed = %g\n'                           %(maxpos, maxspeed))
            logfile1.write('number of particles, approximately = %g\n\n' %(pmax)) 
            logfile1.close() 

            Poincare = np.zeros((pcountmax, 2*pmax))

            p = -1

            for pos in range(1, maxpos): #pos in range(maxpos): current position from 0 to maxpos

                a00 = ( -1.0 + 2.0 * (pos * 1.0 / maxpos) )/ (bn * delta) #chosen x0*delta*bn from -1 exlusive to 1 exclusive
                a01 = 0.0
                a02 = 0.0   

                speed0 = maxspeed/2 - int(( 1.0 - (bn*delta*a00)**2)**0.5 / h)
                if speed0 == 0:
                    speed0 = 1    

                for speed in range(speed0, maxspeed - speed0 + 1):
                    p +=1     

                    a03 = (-1.0 + 2.0 * (speed * 1.0 / maxspeed)) #distribution of Vx                                                                           
                    a04 = - bn * delta * a00
                    
                    if a03**2 + a04**2 >= 1:
                        continue
                    
                    a05 = (1.0 - a03**2 - a04**2)**0.5 #so that the initial energy equals 1.0, nesessarily > 0

                    a0 = np.array([a00, a01, a02, a03, a04, a05])
                           

                    print 'pos =', pos, ', speed =', speed
                    RKt, Poinc, ZeroIndices = RK4_Poinc(f, a0)            

                    #np.savetxt(r'%s\Trajectory_%d.txt' %(folder, p), \
                               #RKt, ['%.4g','%.4g','%.4g','%.4g','%.4g','%.4g','%g',], delimiter=',') 

                    #np.savetxt(r'%s\Poinc_%d.txt' %(folder, p), Poinc, '%.3g', delimiter=',')

                    Poincare[:, 2*p] = Poinc[:, 0]
                    Poincare[:, 2*p + 1] = Poinc[:, 1] 

                    #indicesfile = open(r'%s\ZeroIndices_%d.txt' %(folder, p), 'w')
                    #for i in ZeroIndices: 
                        #indicesfile.write('%s\n' %i)    
                    #indicesfile.close()   


            np.savetxt(r'%s\Poincare.txt' %(folder), Poincare, '%.3g', delimiter=',')
            np.savetxt(r'%s\Poincare_N.txt' %(folder), Poincare[::2], '%.3g', delimiter=',')
            np.savetxt(r'%s\Poincare_S.txt' %(folder), Poincare[1::2], '%.3g', delimiter=',')

            logfile3 = open(r'%s\log.txt' %(folder), 'w')
            logfile3.close()

