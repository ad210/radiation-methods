import numpy as np
from scipy.integrate import odeint
import math

import tools

# Define spacially and temporally dependent E field by components
def E_fun(r,t,params):
    x, y, z = r
    q, m, L, A, c, N = params
    Ex = 0.
    Ey = 0.
    Ez = 0.
    return np.array([Ex,Ey,Ez]) 

# Define spacially and temporally dependent B field by components
def B_fun(r,t,params):
    x, y, z = r
    q, m, L, A, c, N = params
    Bx = 0.
    By = 0.
    Bz = 0.

    k = .5; d = .25;
    A = -A/4; w = 20;
   
    
    a = (k/8); b = ((5.0-k)/8);
    
    x0 = 0*(a+b) + d;
    x1 = 1*(a+b) + d;    
    x2 = 2*(a+b) + d; 
    x3 = 3*(a+b) + d;

    #Bz =      A*( math.exp( w*(x-x0) ) / (math.exp( w*(x-x0) ) + math.exp( -w*(x-x0) )) - math.exp( w*(x-x0-a) ) / (math.exp( w*(x-x0-a) ) + math.exp( -w*(x-x0-a) )))
    #Bz = Bz - A*( math.exp( w*(x-x1) ) / (math.exp( w*(x-x1) ) + math.exp( -w*(x-x1) )) - math.exp( w*(x-x1-a) ) / (math.exp( w*(x-x1-a) ) + math.exp( -w*(x-x1-a) )))
    #Bz = Bz - A*( math.exp( w*(x-x2) ) / (math.exp( w*(x-x2) ) + math.exp( -w*(x-x2) )) - math.exp( w*(x-x2-a) ) / (math.exp( w*(x-x2-a) ) + math.exp( -w*(x-x2-a) )))
    #Bz = Bz + A*( math.exp( w*(x-x3) ) / (math.exp( w*(x-x3) ) + math.exp( -w*(x-x3) )) - math.exp( w*(x-x3-a) ) / (math.exp( w*(x-x3-a) ) + math.exp( -w*(x-x3-a) )))
    
    Bz = 15*A

    #if x > x0 and x <= x0 + a:
    #    Bz = A;
    #if x > x1 and x <= x1 + a:
    #    Bz = -A;
    #if x > x2 and x <= x2 + a:
    #    Bz = -A;
    #if x > x3 and x <= x3 + a:
    #    Bz = A;
    

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
SHOW_PARAM = 1 
SAVE_FILE = 0
num_particles = 1

for n in range(0,num_particles):
    [q_e,m_e,c,tau] = tools.get_constants()
    
    # ====== Define Undulator Properties ====== #
    K_und = 1.038                   # undulator parameter
    lambda_und = .1106              # spactial period
    N_und = 7.                      # number of periods
    L_und = 3.0 #lambda_und*N_und
    lambda_axs = 2.2E-6             # on-axis wavelength
    gamma_e = 95.69                # electron lorentz factor
    B_und = abs(tau*m_e*c*K_und/(q_e*lambda_und)) # b-field calculation
    undulator_params = [K_und, L_und, lambda_und, N_und, lambda_axs, gamma_e, B_und]
    
    # ======= Initial Values ====== #
    vi = math.sqrt(1-(1./gamma_e)**2)*c
    rx, ry, rz =[0., 0. ,0.]
    vx_init, vy_init, vz_init = [vi,0,0]
    vx, vy, vz = [vx_init,vy_init,vz_init]    

    params = [q_e,m_e,lambda_und,B_und,c,N_und] 
    y_0 = [rx, ry, rz, vx, vy, vz]
    
    # ====== Run Integrater ====== #    
    steps = 2000.
    time_intv = L_und/(vi*steps)
    time_undulat = L_und/vi
    time_range = np.arange(0.0,time_undulat,time_intv)
    usoln = odeint(relativistic_f, y_0, time_range, args=(params,),hmax=time_intv)
    
    # ====== Generate Padding ====== #
    p_length = 0*L_und
    time_pad = p_length/(vi)
    pad_step = time_pad/time_intv
    
    psoln = tools.generate_padding([vx_init,vy_init],time_intv,pad_step)  
    for row in usoln:
        psoln = np.append(psoln,[row],axis=0)
    
    
    # ====== Display and File Saving ====== #
    if SHOW_PARAM:
        tools.print_parameters(undulator_params)   
    if SHOW_PLOT:
        tools.plot_trajectory(psoln.transpose())
    if SAVE_FILE:
        filename = "../clara2-dev/src/data/trace_" + str(n).zfill(4) + ".txt"
        tools.write_csv(filename,psoln,time_intv,0)
