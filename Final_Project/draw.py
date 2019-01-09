#!/usr/bin/env python3

from scipy.spatial.distance import pdist, squareform

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation
from map import Configuration, CreatureType, Map

class SimulationBox:
    def __init__(self, bounds = [-1, Configuration.N, -1, Configuration.M], time_steps = 1, config=Configuration):
        self.map = Map(config=config)
        self.bounds = [-1, config.N, -1, config.M]
        self.state = np.zeros((config.N,config.M)).astype(int)
        self.time_steps = time_steps

    def step(self):
        self.map.simulate()
        self.state = self.map.map
    
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

class Simulation:
    def __init__(self, config=Configuration, time_steps=100):
        self.box = SimulationBox(config=config, time_steps=time_steps)
        #------------------------------------------------------------
        # set up figure and animation
        self.fig = plt.figure()
        self.fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        self.ax = self.fig.add_subplot(111, aspect='equal', autoscale_on=False,
                            xlim=(-1, config.N), ylim=(-1, config.M))


        self.predators, = self.ax.plot([], [], 'bo', ms=10, color='r')
        self.preys, = self.ax.plot([], [], 'bo', ms=10, color='b')

        # rect is the box edge
        self.rect = plt.Rectangle(self.box.bounds[::2],
                            self.box.bounds[1] - self.box.bounds[0],
                            self.box.bounds[3] - self.box.bounds[2],
                            ec='none', lw=2, color='black')
        self.ax.add_patch(self.rect)

    def init_anim(self):
        """initialize animation"""
        self.predators.set_data([], [])
        self.preys.set_data([], [])
        self.rect.set_edgecolor('none')
        return self.predators, self.preys, self.rect

    def animate(self, i):
        """perform animation step"""
        self.box.step()
        
        # update pieces of the animation
        self.rect.set_edgecolor('k')
        
        predator_pos, prey_pos = self.box.position_for_drawing()

        self.predators.set_data(predator_pos[0], predator_pos[1])
        self.preys.set_data(prey_pos[0], prey_pos[1])

        return self.predators, self.preys, self.rect

    def start_simulation(self, should_show_animation=False):
        ani = animation.FuncAnimation(self.fig, self.animate, frames=self.box.time_steps,
                              interval=20, blit=True, init_func=self.init_anim)

        ani.save('./simulation.mp4', fps=5.0, dpi=200, writer='ffmpeg')

        if should_show_animation:
            plt.show()