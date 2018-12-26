#!/usr/bin/env python3

import map

creature_map = map.Map(type=map.Map_Type.WAVE)
creature_map.show_map()
creature_map.simulate(1000)
creature_map.show_map()