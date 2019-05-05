import os 

os.chdir("./clara2-dev/src")
os.system("./executable")
os.system("./process_data")
os.system("../tools/plotRadiation ./my_spectrum_all_000.dat")
