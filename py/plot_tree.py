"""
Code to create the Tree Plots such as in the README. Filepaths are hard-coded.

Adjust the settings in USER PARAMS as desired to create different plots, though feel
free to tinker with the rest of the code if you like!
"""

import matplotlib.pyplot as plt
import numpy as np

########################################################################################
# USER PARAMS

# Specify input filepath here (hard-coded)
file_in = "../data/tree_n100_dmax90_seed0_iseed0.csv"

# View the file in interactive mode?
interactive = True

# Save the file to disk?
# N.B. looking at the plot interactively will break the exported plot.
# The script will check to ensure interactive mode is off before saving.
save_file = False
########################################################################################

# Load data
data = np.loadtxt(file_in, delimiter = ",")
x = data[:,0]
y = data[:,1]

# Generate figure
fig, ax = plt.subplots(facecolor = "black")
ax.set_xlim(min(x) - 2, max(x) + 2)
ax.set_ylim(min(y) - 2, max(y) + 2)

# Source: https://stackoverflow.com/questions/48172928/scale-matplotlib-pyplot-axes-scatter-markersize-by-x-scale/48174228#48174228
# Scatter class which will dynamically rescale the points to a diameter specified
# Looks really good in interactive mode. Would recommend giving it a try.
class scatter():
    def __init__(self,x,y,ax,size=1,**kwargs):
        self.n = len(x)
        self.ax = ax
        self.ax.figure.canvas.draw()
        self.size_data=size
        self.size = size
        self.sc = ax.scatter(x,y,s=self.size,**kwargs)
        self._resize()
        self.cid = ax.figure.canvas.mpl_connect('draw_event', self._resize)

    def _resize(self,event=None):
        ppd=72./self.ax.figure.dpi
        trans = self.ax.transData.transform
        s =  ((trans((1,self.size_data))-trans((0,0)))*ppd)[1]
        if s != self.size:
            self.sc.set_sizes(s**2*np.ones(self.n))
            self.size = s
            self._redraw_later()
    
    def _redraw_later(self):
        self.timer = self.ax.figure.canvas.new_timer(interval=10)
        self.timer.single_shot = True
        self.timer.add_callback(lambda : self.ax.figure.canvas.draw_idle())
        self.timer.start()

# Create the scatter plot and fiddle with the aspect ratio and colours etc.
sc = scatter(x,y,ax, size=2, linewidth=0, color="#FFFFFF")
ax.set_aspect(1)
plt.axis("off")
plt.grid("off")

# Show the plot to the user if desired
if interactive == True:
    plt.show()

# Save to disk if desired, and not if we've looked at the plot interactively
# Adjust the DPI to set the detail of the .png. Higher dpi files are considerably larger.
# 1000 is a decent level of detail for <100,000 particles.
elif save_file == True :
    plt.savefig("../fig/" + file_in[8:-4] + ".png", dpi = 1000)