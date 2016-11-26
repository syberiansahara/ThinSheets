import math
from math import tanh, cos

def Field(x, y, z):
      
    bx0 = 1    
    by0 = 1
    bn = 0.10
    
    norm = (bx0**2 + by0**2 + bn**2)**0.5
    bx0 = bx0 / norm
    by0 = by0 / norm
    bn = bn / norm

    field_components = dict()
    
    field_components['bx0'] = bx0
    field_components['by0'] = by0
    field_components['bz'] = bn
        
    field_components['bx'] = bx0 * tanh(z / delta)
    
    
    #'''
    field_components['fieldtype'] = 'constant'
    field_components['by'] = by0  
    #'''
    
    '''
    field_components['fieldtype'] = 'bell'
    if abs(z) < delta:
        field_components['by'] = by0 * cos( (z / delta) * (pi / 2.0))
    else:
        field_components['by'] = 0.0
    '''
    
    '''
    field_components['fieldtype'] = 'antisymm'
    if abs(z) < delta:
        field_components['by'] = - by0 * sin( z * pi / delta)
    else:
        field_components['by'] = 0.0
    '''
  

    return field_components



    