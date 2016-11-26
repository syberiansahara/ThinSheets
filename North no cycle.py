
# coding: utf-8

# # Diploma invariant

# Importing libraries

# In[1]:

import numpy as np
import scipy as sp
from scipy.integrate import quad, simps
import math
from math import atan, sin, cos, acosh, exp, log, cosh, tanh, pi
from numpy.linalg import norm


# In[2]:

def sign(a):
    if a > 0.0:
        return 1.0
    elif a < 0.0:
        return - 1.0
    else:
        return 0.0


# Declaring constants

# In[3]:

delta = 1.0
eps = 1e-7 #for energy conservation in Runge-Kutta method


# In[4]:

theroot = r'C:\Users\syber\Syber Python\Trial'


# In[5]:

bx0 = 1    
by0 = 0.5
bn = 0.10
thebn = bn
#fieldtype = 'constant'
#fieldtype = 'bell'
fieldtype = 'antisymm'
#fieldtype = 'shearless'


root = r'%s\%s\by=%.2fbx\bn=%.2fbx\North' %(theroot, fieldtype, by0, bn)

norm = (bx0**2 + by0**2 + bn**2)**0.5
bx0 = bx0 / norm
by0 = by0 / norm
bn = bn / norm

def Field(x, y, z):
         
    field_components = dict()
    
    field_components['bx'] = bx0 * tanh(z / delta)  
    field_components['bz'] = bn
    
    if fieldtype == 'constant':
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

# In[6]:

def f(b, a):
        
    q = np.empty((6,))
    q[0] = a[3]
    q[1] = a[4]
    q[2] = a[5]
    q[3] = delta * (b['bz']*a[4] - b['by']*a[5])
    q[4] = delta * (b['bx']*a[5] - b['bz']*a[3])
    q[5] = delta * (b['by']*a[3] - b['bx']*a[4]) 
    
    return q


# Function returning invariant for all trajectory points

# In[9]:

def InvarNum(RKt, ZeroVz):
    Inv = np.zeros((RKt.shape[0], 1))
 
    if len(ZeroVz) % 2 == 0:
        del ZeroVz[-1]
        
    for i in range((len(ZeroVz) - 1) / 2):
        i1 = ZeroVz[2*i]        
        i2 = ZeroVz[2*i + 1]     
        i3 = ZeroVz[2*i + 2]    
        
        y1 = RKt[i1:i2, 5]
        x1 = RKt[i1:i2, 2]
        
        integr1 = simps(y1, x1, even='first')
        int1 = integr1
        
        y2 = RKt[i2:i3, 5]
        x2 = RKt[i2:i3, 2]
        
        integr2 = simps(y2, x2, even='first')
        int2 = integr2
        
        Inv[i1:i3, 0] = int1 + int2
        
    Inv[0:ZeroVz[0], 0] = Inv[ZeroVz[0]+2, 0]
    Inv[ZeroVz[-2]:RKt.shape[0], 0] = Inv[ZeroVz[-1]-2, 0]
    return Inv


# Function returning array with invariant surges

# In[1]:

def SurgeFunction(ZeroIndices, ZeroVz, Invariant, RKt):
    result = []
    if len(ZeroIndices) == 0:
        return result
    else:
        zi1 = ZeroIndices[0]        
        if len(ZeroIndices) == 1:
            zi2 = Invariant.shape[0] - 1
        else:                        
            for k in range(1, len(ZeroIndices)):
                        
                Serpentine = True
                
                if len([ind for ind in ZeroVz if ind > ZeroIndices[k-1] and ind <= ZeroIndices[k]]) >= 2:
                    Serpentine = False
                
                if not Serpentine:
                    break
                
            if Serpentine:
                zi2 = ZeroIndices[k]
            elif k == 1:
                zi2 = ZeroIndices[1]
            else:
                zi2 = ZeroIndices[k-1]
            
        mean1 = np.mean(Invariant[0 : zi1]) 
        mean2 = np.mean(Invariant[zi1 : zi2])
        
        result.append((mean2 - mean1) / mean1)
        
        zi1zi2 = []
        zi1zi2.append(zi1)
        zi1zi2.append(zi2)
               
    
        return result, zi1zi2


# In[11]:

def InitialConditions(a02, phase, pitch): #angle setting for Invariant calculations
        
    b = Field(0, 0, a02)

    if b['bx'] == 0 and b['by'] == 0:
        alpha = (pi / 2)
    else:
        alpha = atan(b['bz'] / (b['bx']**2 + b['by']**2)**0.5)  

    if b['bx'] == 0:
        gamma = (pi / 2) * sign(b['by'])
    elif b['bx'] > 0:
        gamma = atan(b['by'] / b['bx']) 
    else: 
        gamma = pi + atan(b['by'] / b['bx'])
        
    phi = (2.0 * pi) * (phase * 1.0) / maxphase
    
    theta = (1.4) * (pitch * 1.0) / maxpitch
    
    
    
    if a02 > 0.0:
        theta = pi - theta #so that the particle moves towards the neutral sheet
    
    print 'phi =', round(phi, 3), ', theta =', round(theta, 3)
    
    a03 = (-sin(theta) * sin(alpha) * cos(phi) + cos(theta) * cos(alpha)) * cos(gamma) -         (sin(theta) * sin(phi)) * sin(gamma)
    a04 = (-sin(theta) * sin(alpha) * cos(phi) + cos(theta) * cos(alpha)) * sin(gamma) +         (sin(theta) * sin(phi)) * cos(gamma)
    a05 = cos(theta) * sin(alpha) + sin(theta) * cos(alpha) * cos(phi)

        
    a0 = np.array([a00, a01, a02, a03, a04, a05])
    
    return a0    


# Defining Runge-Kutta method

# In[12]:

def RK4_Inv(f, a0):
    Cy = a0[4]  + delta*bn*a0[0] - delta**2 * bx0 * log( cosh( a0[2] / delta ))
    a0[0] = a0[0] - Cy / (bn*delta) #shift of coordinates    
    
    
    RK = a0.reshape((1, 6)) #trajectory array
    times = [0.0] #array with times
    Poinc = np.zeros((pcountmax, 2)) #array with found poincare points for the current particle
    ZeroIndices = [] #indices of poincare points in trajectory array
    ZeroVz = [] #indices of null Vz speed
    
    en0 = (a0[3]**2 + a0[4]**2 + a0[5]**2)**0.5 #initial energy equals 1.0
    
    t = 0.0 #initial time
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
        
        
        if OldRK[5] * NewRK[5] < 0: #Vz speed sign change condition
            ZeroVz.append(i)
            
        ##=================== Poincare==================
        
        if OldRK[2] * NewRK[2] < 0: #z coordinate sign change condition
            
            ZeroIndices.append(i)
            
            if pcount < pcountmax:        
                                    
                #linear approximation of bn*delta*(x)
                Poinc[pcount, 0] = bn * delta * (OldRK[0] + OldRK[2] * (NewRK[0] - OldRK[0]) /                                                  (OldRK[2] - NewRK[2]))
                
                #linear approximation of Vx
                Poinc[pcount, 1] = OldRK[3] + OldRK[2] * (NewRK[3] - OldRK[3]) / (OldRK[2] - NewRK[2])
                                    
                pcount +=1
        
        if len(ZeroIndices) >= 2 and ( len([ind for ind in ZeroVz if ind > ZeroIndices[-1]]) >= 2                                       or len([ind for ind in ZeroVz if ind > ZeroIndices[-2] and ind <= ZeroIndices[-1]]) >= 2):
            break #Serpentine mode is over          
        
        if pcount == pcountmax:
            break
        
        ##==============================================

        if abs(NewRK[2] / a02) > 2.0: #particle too far from the layer
            toofar = True
            break
        
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
   
    return RKt, Poinc, ZeroIndices, ZeroVz


# Setting

# In[13]:

dt0 = 0.01
tfinal = 800.0
pcountmax = 19 # maximum Poincare points, should be odd
maxphase = 10 #number of different pitch angles
maxpitch = 11 #number of different phases + 1
pmax = maxphase * (maxpitch - 1)

a00 = 0.0
a01 = 0.0  


# ## One particle

# In[ ]:

a02 = 5.0
    
logfile1 = open(r'%s\parameters.txt' %(root), 'w')
logfile1.write('Field type = %s\nbx0 = %g\nby0 = %g\nbn = %g\n\ndelta = %g\ndt0 = %g\ntfinal = %g\npcountmax = %g\n\n'               %(fieldtype, bx0, by0, bn, delta, dt0, tfinal, pcountmax))
logfile1.write('maxphase = %g\nmaxpitch = %g\n'               %(maxphase, maxpitch))
logfile1.write('number of particles = %g\n\n' %(pmax)) 
logfile1.close()  

phase = 0
pitch = 1

p = phase * (maxpitch - 1) + pitch - 1

a0 = InitialConditions(a02, phase, pitch)

RKt, Poinc, ZeroIndices, ZeroVz = RK4_Inv(f, a0)

InvariantNum = InvarNum(RKt, ZeroVz)
SurgeNum, zi1zi2 = SurgeFunction(ZeroIndices, ZeroVz, InvariantNum, RKt)

np.savetxt(r'%s\Trajectory_%d.txt' %(root, p),                    RKt, ['%.4g','%.4g','%.4g','%.4g','%.4g','%.4g','%g',], delimiter=',')  
np.savetxt(r'%s\Poinc_%d.txt' %(root, p),                    Poinc, '%.3g', delimiter=',')
np.savetxt(r'%s\InvariantNum_%d.txt' %(root, p),                    InvariantNum, '%.4g', delimiter=',')

indicesfile0 = open(r'%s\ZeroVz_%d.txt' %(root, p), 'w')
for i in ZeroVz: 
    indicesfile0.write('%s\n' %i)    
indicesfile0.close()     
    
indicesfile = open(r'%s\ZeroIndices_%d.txt' %(root, p), 'w')
for i in ZeroIndices: 
    indicesfile.write('%s\n' %i)    
indicesfile.close()

indfile = open(r'%s\zi1zi2_%d.txt' %(root, p), 'w')
for i in zi1zi2: 
    indfile.write('%s\n' %i)    
indfile.close()

SurgeNum = np.array(SurgeNum)
print 'Surge: ', SurgeNum[0]

np.savetxt(r'%s\SurgeNum_%d.txt' %(root, p), SurgeNum, '%g', delimiter=',')


# ## Multiple particles (with invariant)

# In[ ]:

a02 = 5

logfile1 = open(r'%s\parameters.txt' %(root), 'w')
logfile1.write('Field type = %s\nbx0 = %g\nby0 = %g\nbn = %g\n\ndelta = %g\ndt0 = %g\ntfinal = %g\npcountmax = %g\n\n'               %(fieldtype, bx0, by0, bn, delta, dt0, tfinal, pcountmax))
logfile1.write('maxphase = %g\nmaxpitch = %g\n'               %(maxphase, maxpitch))
logfile1.write('number of particles = %g\n\n' %(pmax)) 
logfile1.close() 

Poincare = np.zeros((pcountmax, 2*pmax))

AllSurgesNum = []

p = -1

for phase in range(maxphase):
    for pitch in range(1, maxpitch):
        p +=1 
        
        print 'phase =', phase, ', pitch =', pitch
        a0 = InitialConditions(a02, phase, pitch)        
        
        RKt, Poinc, ZeroIndices, ZeroVz = RK4_Inv(f, a0)
        
        InvariantNum = InvarNum(RKt, ZeroVz)
        
        SurgeNum, zi1zi2 = SurgeFunction(ZeroIndices, ZeroVz, InvariantNum, RKt)  

        AllSurgesNum.extend(SurgeNum)
        
        Poincare[:, 2*p] = Poinc[:, 0]
        Poincare[:, 2*p + 1] = Poinc[:, 1] 
           
        #np.savetxt(r'%s\Trajectory_%d.txt' %(root, p), \
                   #RKt, ['%.4g','%.4g','%.4g','%.4g','%.4g','%.4g','%g',], delimiter=',') 
        
        #indicesfile0 = open(r'%s\ZeroVz_%d.txt' %(root, p), 'w')
        #for i in ZeroVz: 
            #indicesfile0.write('%s\n' %i)    
        #indicesfile0.close()   
        
        #indicesfile = open(r'%s\ZeroIndices_%d.txt' %(root, p), 'w')
        #for i in ZeroIndices: 
            #indicesfile.write('%s\n' %i)    
        #indicesfile.close() 
        
        #indfile = open(r'%s\zi1zi2_%d.txt' %(root, p), 'w')
        #for i in zi1zi2: 
            #indfile.write('%s\n' %i)    
        #indfile.close()

        
        #np.savetxt(r'%s\InvariantNum_%d.txt' %(root, p), \
                   #InvariantNum, '%.4g', delimiter=',')
        
        #np.savetxt(r'%s\Poinc_%d.txt' %(root, p), \
                   #Poinc, '%.3g', delimiter=',')
            
            
AllSurgesNum = np.array(AllSurgesNum)

print 'Surges Numerical. Mean: ', np.mean(AllSurgesNum), '. Mean of **2:', np.mean(AllSurgesNum**2)
    
np.savetxt(r'%s\Poincare.txt' %(root), Poincare, '%.3g', delimiter=',')
np.savetxt(r'%s\AllSurgesNum.txt' %(root), AllSurgesNum, '%g', delimiter=',')

logfile3 = open(r'%s\meansnum.txt' %(root), 'w')
logfile3.write('%g,%g,%g,%g,%g' %(thebn, np.mean(AllSurgesNum), np.mean(AllSurgesNum**2), log(thebn), log(np.mean(AllSurgesNum**2))))
logfile3.close()

logfile4 = open(r'%s\log.txt' %(root), 'w')
logfile4.close()


# p = phase * (maxpitch - 1) + pitch - 1
