import os,csv,math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from decimal import Decimal


# === script parameters === #
SAVE_FILE = True
SHOW_PLOT = True

full_array = [];    # Saves spectrum at all trace points 
step_size  = 2      # size of each splice
num_points = 3000   # total number of points in trace file


# === file handling === #
working_dir = "../clara2-dev/src"
os.chdir(working_dir)
os.system("cp data/trace_0000.txt data/temp.txt")

temp_file = "data/temp.txt"
trace_file = "data/trace_0000.txt"
trace_array = np.loadtxt(open(temp_file,'r'),delimiter = " ")



# === main loop === #
for cut_off in range(0, 200, step_size):
    
    print "\n\nINDEX: " + str(cut_off) +" out of " +str(num_points)+"\n"

    temp = np.array(trace_array[:-(num_points-cut_off)])
    np.savetxt(trace_file, temp, delimiter=" ")

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
    full_array.append(summed_freqs.tolist())


# === remove temp file === #
os.system("cp data/temp.txt data/trace_0000.txt")
os.system("rm data/temp.txt")

# === save and display === #
full_array = np.array(full_array).transpose()
if SHOW_PLOT:
    fig = plt.figure("Radiation Intensity (Angle vs. wavelength)")
    ax = fig.gca()
    ax.matshow(full_array)
    #ax.matshow(np.diff(full_array))
    plt.show()

if SAVE_FILE:
    np.savetxt("../../splice-script/radiation_splice.csv", full_array, delimiter=",")


