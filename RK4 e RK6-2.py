import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

def f(X,Y,Z):
    return Z

def g(X,Y,Z):
    return Z**2 - 1

def rk4(a,b,h,y0,z0):
    n = int((b-a)/h)
    x = np.zeros(n+1)
    y = np.zeros(n+1)
    z = np.zeros(n+1)
    x[0] = a
    y[0] = y0
    z[0] = z0
    
    for k in range(0, n):
        k1 = f(x[k], y[k], z[k])*h
        l1 = g(x[k], y[k], z[k])*h
        
        k2 = f(x[k]+h/2, y[k]+k1/2, z[k]+l1/2)*h
        l2 = g(x[k]+h/2, y[k]+k1/2, z[k]+l1/2)*h
        
        k3 = f(x[k]+h/2, y[k]+k2/2, z[k]+l2/2)*h
        l3 = g(x[k]+h/2, y[k]+k2/2, z[k]+l2/2)*h
        
        k4 = f(x[k]+h, y[k]+k3, z[k]+l3)*h
        l4 = g(x[k]+h, y[k]+k3, z[k]+l3)*h
        
        
        y[k+1] = y[k] + (1/6)*(k1+2*k2+2*k3+k4)
        z[k+1] = z[k] + (1/6)*(l1+2*l2+2*l3+l4)
        x[k+1] = x[k] + h
        
    
    xreal = np.arange(a,b,0.001)
    yreal = xreal - np.log(sc.e**(2*xreal) + sc.e**2) + np.log(1 + sc.e**2)
    
    plt.figure(dpi=1000)
    plt.plot(x,y,'--bx')
    plt.plot(xreal,yreal,'r')
    
    print(x)
    print(y)
    #print(z)

def rk6(a,b,h,y0,z0):
    n = int((b-a)/h)
    x = np.zeros(n+1)
    y = np.zeros(n+1)
    z = np.zeros(n+1)
    x[0] = a
    y[0] = y0
    z[0] = z0
    
    for k in range(0, n):
        k1 = f(x[k], y[k], z[k])*h
        l1 = g(x[k], y[k], z[k])*h
        
        k2 = f(x[k]+h/3, y[k]+k1/3, z[k]+l1/3)*h
        l2 = g(x[k]+h/3, y[k]+k1/3, z[k]+l1/3)*h
        
        k3 = f(x[k]+h*2/3, y[k]+k2*2/3, z[k]+l2*2/3)*h
        l3 = g(x[k]+h*2/3, y[k]+k2*2/3, z[k]+l2*2/3)*h
        
        k4 = f(x[k]+h/3, y[k]+(k1/12+k2/3-k3/12), z[k]+(l1/12+l2/3-l3/12))*h
        l4 = g(x[k]+h/3, y[k]+(k1/12+k2/3-k3/12), z[k]+(l1/12+l2/3-l3/12))*h
        
        k5 = f(x[k]+h/2, y[k]+(-k1/10+9*k2/8-3*k3/16-3*k4/8), z[k]+(-l1/10+9*l2/8-3*l3/16-3*l4/8))*h
        l5 = g(x[k]+h/2, y[k]+(-k1/10+9*k2/8-3*k3/16-3*k4/8), z[k]+(-l1/10+9*l2/8-3*l3/16-3*l4/8))*h
        
        k6 = f(x[k]+h/2, y[k]+(9*k2/8-3*k3/8-3*k4/4+k5/2), z[k]+(9*l2/8-3*l3/8-3*l4/4+l5/2))*h
        l6 = g(x[k]+h/2, y[k]+(9*k2/8-3*k3/8-3*k4/4+k5/2), z[k]+(9*l2/8-3*l3/8-3*l4/4+l5/2))*h
        
        y[k+1] = y[k] + (11*k1/120 + 27*k3/40 + 27*k4/40 - 4*k5/15 - 4*k6/15)
        z[k+1] = z[k] + (11*l1/120 + 27*l3/40 + 27*l4/40 - 4*l5/15 - 4*l6/15)
        x[k+1] = x[k] + h
        
    
    xreal = np.arange(a,b,0.001)
    yreal = xreal - np.log(sc.e**(2*xreal) + sc.e**2) + np.log(1 + sc.e**2)
    
    plt.plot(x,y,'--go')
    plt.plot(xreal,yreal,'r')
    
    print(x)
    print(y)
        
z0 = (sc.e**2-1)/(sc.e**2+1)

rk4(0, 10, 1, 0, z0)
rk6(0, 10, 1, 0, z0)


plt.savefig('plot.png')


