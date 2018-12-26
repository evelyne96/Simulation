from enum import IntEnum
import numpy as np
import random

class CreatureType(IntEnum):
    NOTHING = 0
    PREDATOR = 1
    PREY = 2

class Map_Type(IntEnum):
    RANDOM = 0
    BLOB = 1
    WAVE = 2

class Configuration():
    PRAY_DEATH_RATE     = 0.2
    PREDATOR_DEATH_RATE = 0.1
    PREDATOR_BIRTH_RATE = 0.8

    PREY_NR             = 225
    PREDATOR_NR         = 300


class Map():

    def __init__(self, N=30, M=30, type=Map_Type.RANDOM):
        self.N = N
        self.M = M
        self.map = np.zeros((N,M)).astype(int)
        self.initialize_map(type)

    def initialize_map(self,type):
        if type == Map_Type.RANDOM:
            self.init_random()
        elif type == Map_Type.BLOB:
            self.init_blob()
        else:
            self.init_wave()

    def init_random(self):
        self.init_creature_randomly(Configuration.PREDATOR_NR, CreatureType.PREDATOR)
        self.init_creature_randomly(Configuration.PREY_NR, CreatureType.PREY)
    
    def init_creature_randomly(self, number, creatureType):
        for i in range(0, number):
            posX = random.randint(0, self.N-1)
            posY = random.randint(0, self.M-1)
            while (self.map[posX][posY] != CreatureType.NOTHING):
                posX = random.randint(0, self.N-1)
                posY = random.randint(0, self.M-1)
            self.map[posX][posY] = creatureType

    def init_blob(self):
        self.map[10][10:15] = CreatureType.PREY
        for i in range(11,15):
            self.map[i][10:14] = CreatureType.PREY
        for i in range(11,15):
            self.map[i][14] = CreatureType.PREDATOR

    def init_wave(self):
        for i in range(0,self.N):
            self.map[i][13] = CreatureType.PREY
            self.map[i][14] = CreatureType.PREY
            self.map[i][15] = CreatureType.PREDATOR

    def show_map(self):
        print(self.map)