import numpy as np
from scipy.integrate import odeint
import math
import matplotlib.pyplot as plt

def odeint_eq(coords,time_range,params):
    rx, ry, px, py = coords
    B, dt = params
    
    derivs = [px + B*ry/2,py - B*rx/2,B*py/2 - rx*(B**2)/4,-B*px/2 - ry*(B**2)/4]
    return derivs

def euler(coords,time_range,params):
    rx, ry, px, py = coords
    B, dt = params
    qsoln = [[rx],[ry]]
    psoln = [[px],[py]]

    esoln = [(1/2)*(px**2 + py**2) - (B/2)*(py*rx-px*ry) + (B**2/8)*(rx**2 + ry**2)]
    for t in time_range:
        rx_n = rx + (px + B*ry/2)*dt
        ry_n = ry + (py - B*rx/2)*dt
        px_n = px - (-B*py/2 + rx*(B**2)/4)*dt
        py_n = py - (+B*px/2 + ry*(B**2)/4)*dt
        rx = rx_n; ry = ry_n;
        px = px_n; py = py_n
        
        qsoln[0].append(rx)
        qsoln[1].append(ry)
        psoln[0].append(px)
        psoln[1].append(py)
        esoln.append((1/2)*(px**2 + py**2) - (B/2)*(py*rx-px*ry) + (B**2/8)*(rx**2 + ry**2))
    return qsoln,psoln,esoln

def symplectic_pos(coords,time_range,params):
    rx, ry, px, py = coords
    B, h = params
    qsoln = [[rx],[ry]]
    psoln = [[px],[py]]
    esoln = [(1/2)*(px**2 + py**2) - (B/2)*(py*rx-px*ry) + (B**2/8)*(rx**2 + ry**2)]
    
    h = h/2
    for t in time_range:
        #rx_1 = (rx + px*h/2 + B*ry*h/4 + (B*py*h**2)/8)/(1 + ((B*h)**2)/16)
        #ry_1 = (ry + py*h/2 - B*rx*h/4 - (B*px*h**2)/8)/(1 + ((B*h)**2)/16)
        #px_2 = px - (-B*py/2 + rx_1*(B**2)/4)*dt
        #py_2 = py - (+B*px/2 + ry_1*(B**2)/4)*dt
        
        rx_1 = (rx + px*h/2 + B*ry*h/4 + (B*py*h**2)/8)/(1 + ((B*h)**2)/16)
        ry_1 = (ry + py*h/2 - B*rx*h/4 - (B*px*h**2)/8)/(1 + ((B*h)**2)/16)
        px_2 = (py + (B*h/4)*(py - B*rx_1/2) + (h/2)*((B/2)*(py - (h/2)*(B*px/2 + (B**2/2)*(ry_1))) - rx*B**2/4))/(1 + (h*B/4)**2) 
        py_2 = (px - (B*h/4)*(px + B*ry_1/2) - (h/2)*((B/2)*(px + (h/2)*(B*py/2 - (B**2/2)*(rx_1))) + ry*B**2/4))/(1 + (h*B/4)**2)
        
        rx_2 = rx_1 + (h/2)*(px_2 + ry_1*B/2)
        ry_2 = ry_1 + (h/2)*(py_2 + rx_1*B/2)

        rx = rx_2; ry = ry_2;
        px = px_2; py = py_2;
        
        qsoln[0].append(rx)
        qsoln[1].append(ry)
        psoln[0].append(px)
        psoln[1].append(py)
        esoln.append((1/2)*(px**2 + py**2) - (B/2)*(py*rx-px*ry) + (B**2/8)*(rx**2 + ry**2))
    
    return qsoln, psoln, esoln

rx, ry = [0,0]
px, py = [1 ,0]

B = 0.1; dt = 0.1;

coords = [rx,ry,px,py]
params = [B,dt]

time_range = np.arange(0,10000*dt,dt)

qsoln, psoln, esoln = symplectic_pos(coords,time_range,params)
osoln = odeint(odeint_eq, coords, time_range, args=(params,)) 

a = list(map(lambda t: math.pi/2 - t/10, time_range))
x = list(map(lambda t: 10*math.cos(t),a))
y = list(map(lambda t: 10*math.sin(t)-10,a))

time = np.arange(0,len(qsoln))
H_ode = list(map(lambda i: osoln[2][i]**2 + osoln[3][i]**2 - (B/2)*(osoln[2][i]*osoln[1][i] - osoln[3][i]*osoln[0][i]) +(B**2 / 8)*(osoln[0][i]**2 + osoln[1][i]**2),time))
H_sym = list(map(lambda i: psoln[0][i]**2 + psoln[1][i]**2 - (B/2)*(psoln[0][i]*qsoln[1][i] - psoln[1][i]*qsoln[0][i]) +(B**2 / 8)*(qsoln[0][i]**2 + qsoln[1][i]**2),time))

fig = plt.figure()
ax = fig.gca()
ax.plot(qsoln[0])
ax.plot(osoln.transpose().tolist()[0])
#ax.plot(y)
#ax.plot(time,osoln.transpose().tolist()[0],time,qsoln[0])
plt.show()

