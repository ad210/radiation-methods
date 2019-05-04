import numpy as np
from scipy.integrate import odeint
import math

import tools

# Define spacially and temporally dependent E field by components
def E_fun(r,t,params):
    x, y, z = r
    Ex = 0.
    Ey = 0.
    Ez = 0.
    return np.array([Ex,Ey,Ez]) 

# Define spacially and temporally dependent B field by components
def B_fun(r,t,params):
    x, y, z = r
    q, m, L, A, c, N = params
    
    Bx = 0.0; By = 0.0; Bz = 0.0;
    if x>= 0 and x < L/2:
        Bz = -((A/2)*math.cos(x*(2*math.pi/L)) - (A/2))
    elif x>=L/2 and x < L*(N+1/2.0):
        Bz = -A*math.cos(x*(2*math.pi/L))
    elif x>=L*(N+1/2.) and x < L*(N+1):
        Bz =-((A/2)*math.cos(x*(2*math.pi/L)) - (A/2))
    else:
        Bz = 0
    return np.array([Bx,By,Bz])
        

def relativistic_f(y,t,params):
    # Unpacking arguments
    rx, ry, rz, vx, vy, vz = y 
    q, m, L, A, c, N = params
    
    v = math.sqrt(vx**2 + vy**2 + vz**2)
    gamma = 1./math.sqrt(1.-(v**2/c**2))

    # Construct vectors
    r = np.array([rx,ry,rz])
    v = np.array([vx,vy,vz])
    E = E_fun(r,t,params)
    B = B_fun(r,t,params)

    # Lorentz force equation
    dev = (1/gamma)*(q/m)*(E + np.cross(v,B) - (v)*(np.dot(v,E)))

    derivs = [vx, vy, vz, dev[0], dev[1], dev[2]]
    return derivs

# ============================ #
# ====== Work Done Here ====== # 
# ============================ #
SHOW_PLOT = 1
SAVE_FILE = 1

[q_e,m_e,c,tau] = tools.get_constants()
q_e = -q_e

def main():
    # ====== Define Undulator Properties ====== #
    K_und = 1.038                   # undulator parameter
    lambda_und = .1106              # spactial period
    N_und = 7.                      # number of periods
    L_und = lambda_und*N_und
    lambda_axs = 2.2E-6             # on-axis wavelength
    gamma_e = 195.69                # electron lorentz factor
    B_und = abs(tau*m_e*c*K_und/(q_e*lambda_und)) # b-field calculation
    undulator_params = [K_und, L_und, lambda_und, N_und, lambda_axs, gamma_e, B_und]

    # ======= Initial Values ====== #
    vi = math.sqrt(1-(1./gamma_e)**2)*c
    v = math.sqrt(2.99778137e+08**2 +  2.49786702e+06**2) 
    vx_init, vy_init = [vi*2.99778137e+08/v, vi*2.49786702e+06/v]
    vx, vy, vz = [vx_init, vy_init, 0]  
    rx, ry, rz = [0., 0., 0.]

    params = [q_e,m_e,lambda_und,B_und,c,N_und] 
    y_0 = [rx, ry, rz, vx, vy, vz]

    # ====== Run Integrater ====== #    
    steps = 2000.
    time_intv = L_und/(vi*steps)
    time_undulat = 5*L_und/vi
    time_range = np.arange(0.0,time_undulat,time_intv)
    usoln = odeint(relativistic_f, y_0, time_range, args=(params,))

    # ====== Generate Padding ====== #
    p_length = 10*L_und 
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
