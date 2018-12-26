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
    PREY_DEATH_RATE     = 0.2
    PREY_BIRTH_RATE     = 0.8
    PREDATOR_DEATH_RATE = 0.1
    PREDATOR_BIRTH_RATE = 0.8

    PREY_NR             = 225
    PREDATOR_NR         = 300


class Map():

    def __init__(self, N=30, M=30, type=Map_Type.RANDOM):
        self.N = N
        self.M = M
        self.map = np.zeros((N,M)).astype(int)
        self.new_map = self.map
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

    def simulate(self, upper_bound):
        """ upper_bound = number of episodes """
        for k in range(0, upper_bound):
            self.new_map = self.map.copy()
            print("Episode: ", k)
            for i in range(0, self.N):
                for j in range(0, self.M):
                    if self.map[i][j] == CreatureType.PREY:
                        self.check_in_all_directions(i, j, CreatureType.PREY)

                    elif self.map[i][j] == CreatureType.PREDATOR:
                        self.check_in_all_directions(i, j, CreatureType.PREDATOR)

            self.kill_starving_predators()
            self.map = self.new_map.copy()

    def check_in_all_directions(self, i , j, creature):
        """ Check in all directions for a given tile"""

        # Up
        self.breed_with_probability(i-1, j, creature)

        # Down
        if i+1 == self.N:
            self.breed_with_probability(0, j, creature)
        else:
            self.breed_with_probability(i+1, j, creature)

        # Right
        if j+1 == self.M:
            self.breed_with_probability(i, 0, creature)
        else:
            self.breed_with_probability(i, j+1, creature)
        
        # Left
        self.breed_with_probability(i, j-1, creature)

        # Up - Right
        if j+1 == self.M:
            self.breed_with_probability(i-1, 0, creature)
        else:
            self.breed_with_probability(i-1, j+1, creature)
        
        # Up - Left
        self.breed_with_probability(i-1, j-1, creature)

        # Down - Left
        if i+1 == self.N:
            self.breed_with_probability(0, j-1, creature)
        else:
            self.breed_with_probability(i+1, j-1, creature)

        # Down - Right
        if i+1 == self.N:
            i=0
        else:
            i=i+1

        if j+1 == self.M:
            self.breed_with_probability(i, 0, creature)
        else:
            self.breed_with_probability(i, j+1, creature)
        
    def breed_with_probability(self, i, j, creature):

        if creature == CreatureType.PREY:
            if self.map[i][j] == CreatureType.NOTHING and np.random.rand() < Configuration.PREY_BIRTH_RATE:
                self.new_map[i][j] = CreatureType.PREY

        elif creature == CreatureType.PREDATOR:
            if self.map[i][j] == CreatureType.PREY and np.random.rand() < Configuration.PREDATOR_BIRTH_RATE:
                self.new_map[i][j] = CreatureType.PREDATOR
    
    def kill_starving_predators(self):
        for i in range(0, self.N):
            for j in range(0, self.M):
                if self.new_map[i][j] == CreatureType.PREDATOR:
                    if self.check_if_its_dead(i, j) == True:
                        self.new_map[i][j] = CreatureType.NOTHING

    def check_if_its_dead(self, i, j):
        """ Return TRUE if there is not any prey in his neighbourhood, FALSE otherwise """

        # Up
        if self.new_map[i-1][j] == CreatureType.PREY:
            return False

        # Down
        if i+1 == self.N:
            if self.new_map[0][j] == CreatureType.PREY:
                return False
        else:
            if self.new_map[i+1][j] == CreatureType.PREY:
                return False

        # Right
        if j+1 == self.M:
            if self.new_map[i][0] == CreatureType.PREY:
                return False
        else:
            if self.new_map[i][j+1] == CreatureType.PREY:
                return False
        
        # Left
        if self.new_map[i][j-1] == CreatureType.PREY:
                return False

        # Up - Right
        if j+1 == self.M:
            if self.new_map[i-1][0] == CreatureType.PREY:
                return False
        else:
            if self.new_map[i-1][j+1] == CreatureType.PREY:
                return False
        
        # Up - Left
        if self.new_map[i-1][j-1] == CreatureType.PREY:
                return False

        # Down - Left
        if i+1 == self.N:
            if self.new_map[0][j-1] == CreatureType.PREY:
                return False
        else:
            if self.new_map[i+1][j-1] == CreatureType.PREY:
                return False

        # Down - Right
        if i+1 == self.N:
            i=0
        else:
            i=i+1

        if j+1 == self.M:
            if self.new_map[i][0] == CreatureType.PREY:
                return False
        else:
            if self.new_map[i][j+1] == CreatureType.PREY:
                return False

        return True
