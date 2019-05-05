import os,csv,math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from decimal import Decimal

working_dir = "../clara2-dev/src"
os.chdir(working_dir)

filename = "data/full_trace.txt"
new_file = "data/trace_0000.txt"
os.system("cp " + new_file + " " + filename)


rad_all_points = [];
dx = 100;
count = 0

for cut_off in range(0,2000,dx):
    print "\n\nRADIATION UP TO " + str(cut_off) +"\n"
    data = open(filename,'r')
    output = open(new_file,'w')

    count = 0;
    for line in data:
        x = float(line.split(' ')[0])
        if x < cut_off*(0.8/2000): #0.002 = 4/2000 (4m chicane)
            output.write(line)

    os.system('./executable --silent')
    os.system('./process_data')


    csv_file= open('my_spectrum_all_000.dat','r')
    csvreader = csv.reader(csv_file,delimiter = '\t')

    spectrum = []
    for row in csvreader:
        row = row[0:len(row)-1]
        row = list(map(lambda x: float(x), row ))
        spectrum.append(list(row))

    summed_freqs = np.array(spectrum).sum(1)
    
    rad_all_points.append(summed_freqs.tolist())

rad_all_points = np.array(rad_all_points).transpose()
fig = plt.figure("Radiation Intensity (Angle vs. wavelength)")
ax = fig.gca()
ax.matshow(np.diff(rad_all_points))
plt.show()

np.savetxt("../../splice_script/radiation_splice.csv", rad_all_points, delimiter=",")
