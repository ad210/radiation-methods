import numpy as np
from scipy.integrate import odeint
import math

import tools

# Define spacially and temporally dependent E field by components
def E_fun(r,v,t,params):
    x, y, z = r
    q, m, c = params
    
    beta = [v[0]/c,v[1]/c,v[2]/c]
    v = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    gamma = 1./math.sqrt(1.-(v**2/c**2))
    
    if x < .5:
        A = -1E-4
    elif x > 5.5 and x < 6.0:
        A = 1E-4
    elif x > 7.0 and x < 7.5:
        A = 1E-4
    elif x > 12.5 and x < 13:
        A = -1E-4
    else:
        A = 0
    
    Ex = 0.
    Ey = gamma*A - (gamma**2)/(gamma+1)*beta[1]*beta[1]*A
    Ez = 0.

    

    return np.array([Ex,Ey,Ez]) 

# Define spacially and temporally dependent B field by components
def B_fun(r,v,t,params):
    x, y, z = r
    q, m, c = params
    
    beta = [v[0]/c,v[1]/c,v[2]/c]
    v = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    gamma = 1./math.sqrt(1.-(v**2/c**2))

    Bx = 0.
    By = 0.
    Bz = 0.
    if x < .5:
        A = -1E-4
    elif x > 5.5 and x < 6.0:
        A = 1E-4
    elif x > 7.0 and x < 7.5:
        A = 1E-4
    elif x > 12.5 and x < 13:
        A = -1E-4
    else:
        A = 0
    return np.array( -gamma*np.cross(beta,[0,A,0]))
        

def relativistic_f(y,t,params):
    # Unpacking arguments
    rx, ry, rz, vx, vy, vz = y 
    q, m, c = params
    
    v = math.sqrt(vx**2 + vy**2 + vz**2)
    gamma = 1./math.sqrt(1.-(v**2/c**2))

    # Construct vectors
    r = np.array([rx,ry,rz])
    v = np.array([vx,vy,vz])
    E = E_fun(r,v,t,params)
    B = B_fun(r,v,t,params)

    # Lorentz force equation
    #dev = (1/gamma)*(q/m)*(E + np.cross(v,B) - (v)*(np.dot(v,E)))
    dev = (1/gamma)*(q/m)*(E + np.cross(v,B))
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
    

    L_und = 1;
    params = [q_e,m_e,c] 
    y_0 = [rx, ry, rz, vx, vy, vz]

    # ====== Run Integrater ====== #    
    steps = 5000.
    time_intv = 10E-12
    time_undulat = steps*time_intv 
    time_range = np.arange(0.0,time_undulat,time_intv)
    usoln = odeint(relativistic_f, y_0, time_range, args=(params,),hmax=time_intv)

    # ====== Generate Padding ====== #
    p_length = 10*steps
    time_pad = p_length*time_intv
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
