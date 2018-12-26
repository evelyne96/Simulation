from enum import Enum
import numpy as np

class CreatureType(Enum):
    PREDATOR = 0
    PREY = 1
    NOTHING = 2

class Creature():

    def __init__(self, creature_type):
        self.type = creature_type

