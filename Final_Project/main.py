#!/usr/bin/env python3

import map
import draw

creature_map = map.Map(type=map.Map_Type.BLOB, filename="blob.csv")
# creature_map.show_map()
creature_map.simulate(100000)
# creature_map.show_map()

# simulation = draw.Simulation()
# simulation.start_simulation()