import math
import numpy as np
import csv 

# ====== Constants ====== #
c = 2.99792458E8
gamma = 350.0
stopping_time = 1E-15 
velocity_init = math.sqrt((1 - math.sqrt(1/gamma))*(c**2))
velocity_stop = 0
stopping_accl = -(velocity_init-velocity_stop)/stopping_time

N_steps = 1000;
time_intv = stopping_time/N_steps;

time = np.arange(time_intv,stopping_time,time_intv)

pos = 0
vel = velocity_init*0.99
acc = stopping_accl

coords = np.matrix([0,0,0,velocity_init,0,0])
for t in time:
    #vel = vel + acc*time_intv
    pos = pos + vel*time_intv
    coords = np.append(coords,np.matrix([pos,0,0,vel,0,0]),axis=0)

# ====== Generate Padding ====== #
p_length = 20*pos
time_pad = p_length/(velocity_init)
pad_step = time_pad/time_intv

x = -p_length
psoln = np.matrix([x,0,0,velocity_init,0,0])
for i in range(0,int(pad_step)):
    x = x + (time_intv*velocity_init)
    psoln = np.append(psoln,[[x,0,0,velocity_init,0,0]],axis=0)

for row in coords:
    psoln = np.append(psoln,row,axis=0)

filename = "trace_" + str(0).zfill(4)
with open("../clara2-dev/src/data/"+filename+".txt","wb") as csvfile:
    writer = csv.writer(csvfile,delimiter = ' ')
    t = 0
    for l in psoln:
        x,y,z,vx,vy,vz = l.tolist()[0]
        x_o, y_o, z_o = [0,0,0];
        r = math.sqrt((x-x_o)**2 + (y-y_o)**2 + (z-z_o)**2)
        bx = (vx/c); by = (vy/c); bz = (vz/c);
        writer.writerow([x,y,z,bx,by,bz,t])
        t = t + time_intv
