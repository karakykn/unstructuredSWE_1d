from source import *

'''Process''' #-1 in inlet means it does not affect the solution
# caseName = "release_in_tank"
# solver = Swe1D(caseName, startTime=0, n=0, inlet=[-1,-1])
# solver.iterativeSolve(CFL=.5, simTime=50, print_step=20)

# caseName = "1dDambreakBenchmark"
# solver = Swe1D(caseName, startTime=0, n=0, inlet=[10,0])
# solver.iterativeSolve(CFL=.5, simTime=50, print_step=20)

caseName = "valveClose"
'''does not work properly, think about the boundaries. also make sure it is the same equation with SWE'''
solver = Swe1D(caseName, startTime=0, n=0, inlet=[-1,-1])
solver.iterativeSolve(CFL=.1, simTime=200, print_step=20)

'''Post-Process'''
solver.plot()