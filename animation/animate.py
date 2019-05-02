from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(5, 8))
my_data = np.genfromtxt('undulator_array.csv', delimiter=',')

#arr = my_data[0:99,:]
#ax.imshow(arr)

def update(i):
    #im_normed = np.random.random((64, 64))
    arr = my_data[(100*i):(100*i+99),:]
    ax.imshow(np.flip(np.transpose(arr),0),vmin=0,vmax=1.5e-30)
    ax.set_axis_off()

anim = FuncAnimation(fig, update, frames=np.arange(0, 199), interval=50)
anim.save('undulator.gif', dpi=80, writer='imagemagick')
plt.close()

