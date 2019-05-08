import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
from mpl_toolkits import mplot3d as p3

num_bodies = 4
time_steps = 5000
time_steps /= 100

def read_data():
	X = []
	Y = []
	Z = []

	with open("output_trajectory.txt","r") as data:
		data = data.read().split('\n')

		for i in range(len(data)-1):
			# traj.append(data[i].split())
			a = data[i].split()
			X.append(a[0])
			Y.append(a[1])
			Z.append(a[2])

	return np.asarray(X, dtype='float'), np.asarray(Y, dtype='float'), np.asarray(Z, dtype='float')

X,Y,Z = read_data()
print("Size of trajectory = "+str(len(X)))

# [width,length, depth]
bounds = [0,199,0,99,0,399] # bounds is size of box : [xmin,xmax,ymin,ymax,zmin,zmax]

# set up figure and animation
fig = plt.figure()
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(bounds[0], bounds[1]), ylim=(bounds[2], bounds[3]))

# particles holds the locations of the particles
particles, = ax.plot([], [], 'bo', ms=6)

# rect is the box edge
rect = plt.Rectangle(bounds[::2],
                     bounds[1] - bounds[0],
                     bounds[3] - bounds[2],
                     ec='none', lw=2, fc='none')
ax.add_patch(rect)

def init():
    """initialize animation"""
    global X,Y, rect
    particles.set_data([], [])
    rect.set_edgecolor('none')
    return particles, rect

def animate(i):
    """perform animation step"""
    global X,Y, rect, ax, fig

    ms = int(fig.dpi * 2 * 0.5 * fig.get_figwidth()
             / np.diff(ax.get_xbound())[0])
    
    # update pieces of the animation
    rect.set_edgecolor('k')
    particles.set_data(X[(i*num_bodies):((i+1)*num_bodies)], Y[(i*num_bodies):((i+1)*num_bodies)])
    particles.set_markersize(ms)
    return particles, rect

ani = animation.FuncAnimation(fig, animate, frames=int(len(X)/num_bodies),
                              interval=200, blit=True, init_func=init)

# plt.show()

ani.save('line.gif', dpi=80, writer='imagemagick')