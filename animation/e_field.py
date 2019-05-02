import numpy as np
import math, csv

# === Linear algebra functions === #
def v_sum(u,v): # vector summation u+v
    ret = list(filter(lambda i: u[i]+v[i] , range(0,len(u))))
    return ret

def v_diff(u,v): # vector subtraction u-v
    ret = list(filter(lambda i: u[i]-v[i] , range(0,len(u))))
    return ret

def s_mult(a,v):
    ret = list(filter(lambda x: a*x, v))
    return ret

# === velocity and radiation field calculations === #
def field_at_point(charge_pos, charge_vel, charge_acc, point_pos, params):
    # Unpacking variables
    [x_0, y_0, z_0] = charge_pos
    [v_x, v_y, z_y] = charge_vel # Charge Velocity normalized to c
    [a_x, a_y, a_z] = charge_acc

    [x_1, y_1, z_0] = point_pos
    [gamma,q,ep_0,c] = params

    # Relavent vectors
    r_ret = math.sqrt( (x_1-x_0)**2 + (y_1-y_0)**2 + (z_1-z_0)**2 )
    
    beta = charge_vel
    beta_dot = charge_acc
    
    beta_mag = np.linalg.norm(beta)
    n = [v_x/beta_mag, v_y/beta_mag, v_z/beta_mag] # direction of travel

    # Equations for near and far-field radiation defined in [2.29] 
    E_vel = s_mult(q/(4*math.pi*ep_0 * r_ret**2 * gamma**2 * (1 - np.dot(beta,n))**3)  ,v_diff(n,beta)) # near field radiation
    E_rad = s_mult(q/(4*math.pi*ep_0 * r_ret * c * (1-np.dot(beta,n))**3),np.cross(n,np.cross(v_diff(n,beta),beta_dot))) # far field radiation

    E_tot = v_sum(E_vel,E_rad)

    return E_tot
    
def point_to_grid(x,y,I,grid,grid_params): # continuous point --> discrete grid
    [dx,dy,num,x_min,y_min] = grid_params

    x_1 = math.floor((x-x_min)/dx); x_2 = math.ceil((x-x_min)/dx)
    y_1 = math.floor((y-y_min)/dy); y_2 = math.ceil((y-y_min)/dy)
    
    cx1 = (x-x_1)/(x_2-x_1); cx2 = 1/cx1;
    cy1 = (y-y_1)/(y_2-y_1); cy2 = 1/cy1;

    grid(x_1,y_1) = grid(x_1,y_1) + cx1*cy1*I
    grid(x_1,y_2) = grid(x_1,y_2) + cx1*cy2*I
    grid(x_2,y_1) = grid(x_2,y_1) + cx2*cy1*I
    grid(x_2,y_2) = grid(x_2,y_2) + cx2*cy2*I
          
    return grid

def main():
    # Some constants
    c = 2.99792E8
    ep_0 = 8.854E-12
    q_e = -1.602E-19
    gamma = 195.69 # for current undulator

        


if __name__ == "__main__":
        main()




