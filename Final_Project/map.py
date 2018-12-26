from enum import Enum
import numpy as np

class CreatureType(Enum):
    PREDATOR = 0
    PREY = 1
    NOTHING = 2

class Map_Type(Enum):
    RANDOM = 0
    BLOB = 1
    WAVE = 2

class Configuration():
    PRAY_DEATH_RATE     = 0.2
    PREDATOR_DEATH_RATE = 0.1
    PREDATOR_BIRTH_RATE = 0.8

class Map():

    def __init__(self, N=30, M=30, type=Map_Type.BLOB):
        self.map = np.zeros((N,M)).astype(int)
        initialize_map(self,type)

    def initialize_map(self,type):
        if type == Map_Type.RANDOM:
            init_random(self)
        elif type == Map_Type.BLOB:
            init_blob(self)
        else:
            init_wave(self)

    def init_random(self):
        pass

    def init_blob(self):
        pass

    def init_wave(self):
        pass

    def show_map(self):
        print(self.map)