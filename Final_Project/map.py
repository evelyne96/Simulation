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
    N                   = 30
    M                   = 30
    PREY_BIRTH_RATE     = 0.8
    PRAY_DEATH_RATE     = 0.2
    PREDATOR_DEATH_RATE = 0.1
    PREDATOR_BIRTH_RATE = 0.8

    PREY_NR             = 225
    PREDATOR_NR         = 300
    MAP_TYPE            = Map_Type.RANDOM


class Map():

    def __init__(self, config=Configuration, filename="statistics.csv"):
        self.config = config
        self.map = np.zeros((self.config.N,self.config.M)).astype(int)
        self.new_map = self.map
        self.initialize_map(config.MAP_TYPE)
        self.statistics_file = filename

    def initialize_map(self,type):
        if type == Map_Type.RANDOM:
            self.init_random()
        elif type == Map_Type.BLOB:
            self.init_blob()
        else:
            self.init_wave()

    def init_random(self):
        self.init_creature_randomly(self.config.PREDATOR_NR, CreatureType.PREDATOR)
        self.init_creature_randomly(self.config.PREY_NR, CreatureType.PREY)
    
    def init_creature_randomly(self, number, creatureType):
        for i in range(0, number):
            posX = random.randint(0, self.config.N-1)
            posY = random.randint(0, self.config.M-1)
            while (self.map[posX][posY] != CreatureType.NOTHING):
                posX = random.randint(0, self.config.N-1)
                posY = random.randint(0, self.config.M-1)
            self.map[posX][posY] = creatureType

    def init_blob(self):
        self.map[10][10:15] = CreatureType.PREY
        for i in range(11,15):
            self.map[i][10:14] = CreatureType.PREY
        for i in range(11,15):
            self.map[i][14] = CreatureType.PREDATOR

    def init_wave(self):
        for i in range(0,self.config.N):
            self.map[i][13] = CreatureType.PREY
            self.map[i][14] = CreatureType.PREY
            self.map[i][15] = CreatureType.PREDATOR

    def show_map(self):
        print(self.map)

    def simulate(self, upper_bound=1):
        """ upper_bound = number of episodes """
        open(self.statistics_file, "w").close()
        for k in range(0, upper_bound):
            self.new_map = self.map.copy()
            num_prey = 0
            num_pred = 0
            for i in range(0, self.config.N):
                for j in range(0, self.config.M):
                    if self.map[i][j] == CreatureType.PREY:
                        num_prey += 1
                        self.check_in_all_directions(i, j, CreatureType.PREY)

                    elif self.map[i][j] == CreatureType.PREDATOR:
                        num_pred += 1
                        self.check_in_all_directions(i, j, CreatureType.PREDATOR)
            with open(self.statistics_file, 'a') as out:
                out.write(str(k) + " " + str(num_pred) + " " + str(num_prey) + '\n')

            self.kill_starving_predators()
            self.map = self.new_map.copy()

    def check_in_all_directions(self, i , j, creature):
        """ Check in all directions for a given tile"""

        # Up
        self.breed_with_probability(i-1, j, creature)

        # Down
        if i+1 == self.config.N:
            self.breed_with_probability(0, j, creature)
        else:
            self.breed_with_probability(i+1, j, creature)

        # Right
        if j+1 == self.config.M:
            self.breed_with_probability(i, 0, creature)
        else:
            self.breed_with_probability(i, j+1, creature)
        
        # Left
        self.breed_with_probability(i, j-1, creature)

        # Up - Right
        if j+1 == self.config.M:
            self.breed_with_probability(i-1, 0, creature)
        else:
            self.breed_with_probability(i-1, j+1, creature)
        
        # Up - Left
        self.breed_with_probability(i-1, j-1, creature)

        # Down - Left
        if i+1 == self.config.N:
            self.breed_with_probability(0, j-1, creature)
        else:
            self.breed_with_probability(i+1, j-1, creature)

        # Down - Right
        if i+1 == self.config.N:
            i=0
        else:
            i=i+1

        if j+1 == self.config.M:
            self.breed_with_probability(i, 0, creature)
        else:
            self.breed_with_probability(i, j+1, creature)
        
    def breed_with_probability(self, i, j, creature):

        if creature == CreatureType.PREY:
            if self.map[i][j] == CreatureType.NOTHING and np.random.rand() < self.config.PREY_BIRTH_RATE:
                self.new_map[i][j] = CreatureType.PREY

        elif creature == CreatureType.PREDATOR:
            if self.map[i][j] == CreatureType.PREY and np.random.rand() < self.config.PREDATOR_BIRTH_RATE:
                self.new_map[i][j] = CreatureType.PREDATOR
    
    def kill_starving_predators(self):
        for i in range(0, self.config.N):
            for j in range(0, self.config.M):
                if self.new_map[i][j] == CreatureType.PREDATOR:
                    if self.check_if_its_dead(i, j) == True:
                        self.new_map[i][j] = CreatureType.NOTHING

    def check_if_its_dead(self, i, j):
        """ Return TRUE if there is not any prey in his neighbourhood, FALSE otherwise """

        # Up
        if self.new_map[i-1][j] == CreatureType.PREY:
            return False

        # Down
        if i+1 == self.config.N:
            if self.new_map[0][j] == CreatureType.PREY:
                return False
        else:
            if self.new_map[i+1][j] == CreatureType.PREY:
                return False

        # Right
        if j+1 == self.config.M:
            if self.new_map[i][0] == CreatureType.PREY:
                return False
        else:
            if self.new_map[i][j+1] == CreatureType.PREY:
                return False
        
        # Left
        if self.new_map[i][j-1] == CreatureType.PREY:
                return False

        # Up - Right
        if j+1 == self.config.M:
            if self.new_map[i-1][0] == CreatureType.PREY:
                return False
        else:
            if self.new_map[i-1][j+1] == CreatureType.PREY:
                return False
        
        # Up - Left
        if self.new_map[i-1][j-1] == CreatureType.PREY:
                return False

        # Down - Left
        if i+1 == self.config.N:
            if self.new_map[0][j-1] == CreatureType.PREY:
                return False
        else:
            if self.new_map[i+1][j-1] == CreatureType.PREY:
                return False

        # Down - Right
        if i+1 == self.config.N:
            i=0
        else:
            i=i+1

        if j+1 == self.config.M:
            if self.new_map[i][0] == CreatureType.PREY:
                return False
        else:
            if self.new_map[i][j+1] == CreatureType.PREY:
                return False

        return True
