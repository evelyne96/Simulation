#!/usr/bin/env python3

import map
import draw
import configparser

# creature_map = map.Map(filename="blob.csv")
# creature_map.show_map()
# creature_map.simulate(100000)
# creature_map.show_map()

config = configparser.ConfigParser()
config.read('config.conf')

sim_conf = map.Configuration()
sim_conf.N = int(config['DEFAULT']['N'])
sim_conf.M = int(config['DEFAULT']['M'])
sim_conf.MAP_TYPE = int(config['DEFAULT']['MAP_TYPE'])

sim_conf.PREY_BIRTH_RATE = float(config['DEFAULT']['PREY_BIRTH_RATE'])
sim_conf.PREDATOR_BIRTH_RATE = float(config['DEFAULT']['PREDATOR_BIRTH_RATE'])
sim_conf.PREDATOR_DEATH_RATE = float(config['DEFAULT']['PREDATOR_DEATH_RATE'])

sim_conf.PREDATOR_NR = int(config['DEFAULT']['PREDATOR_NR'])
sim_conf.PREY_NR = int(config['DEFAULT']['PREY_NR'])

simulation = draw.Simulation(config=sim_conf)
simulation.start_simulation(should_show_animation=True)