import csv, math
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# ====== Universal Constants ====== #
q_e = -1.6021766208E-19     # fundamental charge
m_e = 9.109383561E-31      # mass of electron
c = 2.99792458E8           # speed of light
tau = 2*math.pi            # 2 * pi

def get_constants():
    return [q_e,m_e,c,tau]

def run_odeint(y_0,params,time_intv):
    lambda_und, B_und, N_und = params
    vi = math.sqrt(y_0[0]**2 + y_0[1]**2 + y_0[2]**2)
    time = L_und/vi 
    t = np.arange(0.0,time,time_inv)

def generate_padding(coords,time_intv,pad_step):
    vx, vy = coords
    rx = -vx*time_intv*pad_step
    ry = -vy*time_intv*pad_step

    psoln = np.matrix([rx,ry,0,vx,vy,0])
    for i in range(0,int(pad_step)):
        rx = rx + vx*(time_intv)
        ry = ry + vy*(time_intv)
        pos = np.matrix([rx,ry,0,vx,vy,0])
        psoln = np.append(psoln,pos,axis=0)
    return psoln

def write_csv(filename,isoln,time_intv,t):
    with open(filename,"wb") as csvfile:
        writer = csv.writer(csvfile,delimiter = ' ')
        for l in isoln:
            x,y,z,vx,vy,vz = l.tolist()[0]
            writer.writerow([x,y,z,vx/c,vy/c,vz/c,t])
            t = t + time_intv

def plot_trajectory(psoln,proj="3D"):
    fig = plt.figure()
    if proj == "3D":
        ax = fig.gca(projection='3d')
        ax.plot(psoln[0,:].tolist()[0], psoln[1,:].tolist()[0], psoln[2,:].tolist()[0],'black',label = 'electron path')
    elif proj == "2D":
        ax = fig.gca()
        ax.plot(psoln[0,:].tolist()[0], psoln[1,:].tolist()[0])
    ax.legend()
    plt.show()

def print_parameters(param):
    K_und, L_und, lambda_und, N_und, lambda_axs, gamma_e, B_und = param
    print '====================================='
    print 'Undulator Parameter     : ' + str(K_und)
    print 'Magnetic Field Strength : ' + str(B_und)
    print 'Electron Lorentz Factor : ' + str(gamma_e)
    print 'Spactial Period         : ' + str(lambda_und)
    print 'Length                  : ' + str(L_und)
    print 'Number of Periods       : ' + str(N_und) + "\n"



