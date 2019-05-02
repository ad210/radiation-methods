import os
import numpy as np
import matplotlib.pyplot as plt

filename = 'full_trace.txt'
new_file = '../clara2-dev/src/data/trace_0000.txt'

data = open(filename,'r')
output = open(new_file,'w')

count = 0;
for line in data:
    x = float(line.split(' ')[0])
    if x < 4:
        output.write(line)

os.system('../clara2-dev/src/executable')
os.system('../clara2-dev/src/process_data')
os.system('../clara2-dev/tools/plotRadiation ./my_spectrum_all_000.dat')
