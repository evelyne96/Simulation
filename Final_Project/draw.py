#!/usr/bin/env python3

from scipy.spatial.distance import pdist, squareform

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation
from map import Configuration, CreatureType, Map

class SimulationBox:

    def __init__(self, bounds = [-1, Configuration.N, -1, Configuration.M], time_steps = 1):
        self.bounds = bounds
        self.state = np.zeros((Configuration.N,Configuration.M)).astype(int)
        self.time_steps = time_steps

    def step(self, time_step):
        self.state = Map().map
        # we can get the data for the right timestep here
    
    def position_for_drawing(self):
        predators = [[], []]
        preys = [[], []]
        for posX in range(0, self.state.shape[0]):
            for posY in range(0, self.state.shape[1]):
                if self.state[posX][posY] == CreatureType.PREDATOR:
                    predators[0].append(posX)
                    predators[1].append(posY)
                elif self.state[posX][posY] == CreatureType.PREY:
                    preys[0].append(posX)
                    preys[1].append(posY)
        return predators, preys



dt = 1. / 30 # 30fps

box = SimulationBox()

#------------------------------------------------------------
# set up figure and animation
fig = plt.figure()
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax = fig.add_subplot(111, aspect='equal', autoscale_on=True,
                     xlim=(-1, Configuration.N), ylim=(-1, Configuration.M))


predators, = ax.plot([], [], 'bo', ms=10, color='r')
preys, = ax.plot([], [], 'bo', ms=10, color='b')

# rect is the box edge
rect = plt.Rectangle(box.bounds[::2],
                     box.bounds[1] - box.bounds[0],
                     box.bounds[3] - box.bounds[2],
                     ec='none', lw=2, color='black')
ax.add_patch(rect)


def init():
    """initialize animation"""
    global box, rect
    predators.set_data([], [])
    preys.set_data([], [])
    rect.set_edgecolor('none')
    return predators, preys, rect

def animate(i):
    """perform animation step"""
    global box, rect, dt, ax, fig
    box.step(dt)
    
    # update pieces of the animation
    rect.set_edgecolor('k')
    
    predator_pos, prey_pos = box.position_for_drawing()

    predators.set_data(predator_pos[0], predator_pos[1])
    preys.set_data(prey_pos[0], prey_pos[1])

    return predators, preys, rect

ani = animation.FuncAnimation(fig, animate, frames=600,
                              interval=20, blit=False, init_func=init)

plt.show()