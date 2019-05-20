import numpy as np
from scipy.integrate import odeint
import math

import tools

# Define spacially and temporally dependent E field by components
def E_fun(r,t,params):
    x, y, z = r
    q, m, c = params
    Ex = 0.
    Ey = 0.
    Ez = 0.
    

    d = 0.5
    A = 1E6 #4.576E5

    if x < d:
        Ey = -A
        #Ey =  9.262E6 
    if x > 5.0 and x < 5.0 + d:
        Ey = A 
    if x > 6.0 and x < 6.0 + d:
        Ey = A
    if x >11.0 and x < 11.0+ d:
        Ey = -A

    return np.array([Ex,Ey,Ez]) 

# Define spacially and temporally dependent B field by components
def B_fun(r,t,params):
    x, y, z = r
    q, m, c = params
    Bx = 0.
    By = 0.
    Bz = 0.

    return np.array([Bx,By,Bz])
        

def relativistic_f(y,t,params):
    # Unpacking arguments
    rx, ry, rz, vx, vy, vz = y 
    q, m, c = params
    
    v = math.sqrt(vx**2 + vy**2 + vz**2)
    
    gamma = 1./math.sqrt(1.-(v/c)**2)
    
    # Construct vectors
    r = np.array([rx,ry,rz])
    v = np.array([vx,vy,vz])
    E = E_fun(r,t,params)
    B = B_fun(r,t,params)

    # Lorentz force equation
    #dev = (1/gamma)*(q/m)*(E + np.cross(v,B) - (v)*(np.dot(v,E)))
    dev = (1/gamma)*(q/m)*E
    derivs = [vx, vy, vz, dev[0], dev[1], dev[2]]
    return derivs

SHOW_PLOT = 1
SAVE_FILE = 1

def main():
    [q_e,m_e,c,tau] = tools.get_constants()

    # ======= Initial Values ====== #
    gamma_e = 195.69                    # 100MeV electrons
    vi = math.sqrt(1-(1./gamma_e)**2)*c
    rx, ry, rz =[0., 0. ,0.]
    vx_init, vy_init, vz_init = [vi,0,0]
    vx, vy, vz = [vx_init,vy_init,vz_init]    

    L = 12
    params = [q_e,m_e,c] 
    y_0 = [rx, ry, rz, vx, vy, vz]

    # ====== Run Integrater ====== #    
    steps = 2000.
    time_intv = L/(vi*steps)
    time_undulat = L/vi
    time_range = np.arange(0.0,time_undulat,time_intv)
    usoln = odeint(relativistic_f, y_0, time_range, args=(params,),hmax=time_intv)

    # ====== Generate Padding ====== #
    p_length = 10*L
    time_pad = p_length/(vi)
    pad_step = time_pad/time_intv

    psoln = tools.generate_padding([vx_init,vy_init],time_intv,pad_step)  
    for row in usoln:
        psoln = np.append(psoln,[row],axis=0)

    # ====== Display and File Saving ====== #
    if SHOW_PLOT:
        tools.plot_trajectory(psoln.transpose())
    if SAVE_FILE:
        filename = "../clara2-dev/src/data/trace_" + str(0).zfill(4) + ".txt"
        tools.write_csv(filename,psoln,time_intv,0)

if __name__ == "__main__":
    main()
