import cProfile, pstats, sys, os
from ipi.engine.simulation import Simulation
sim = Simulation.load_from_xml(open("input.xml"))
pr = cProfile.Profile(); pr.enable()
try: sim.run()
finally:
    pr.disable(); pr.dump_stats("prof.out")
os._exit(0)
