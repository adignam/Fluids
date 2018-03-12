# Fluids 1D Sod shock

from riemannsolver import RiemannSolver
import matplotlib.pyplot as plt


# The Cell class
# This sets up the parameters of the cell
class Cell(object):
    def __init__(self):
        self._volume = 0.
        self._mass = 0.
        self._momentum = 0.
        self._energy = 0.
        self._density = 0.
        self._velocity = 0.
        self._pressure = 0.
        self._right_ngb = None
        self._surface_area = 1.
        
# The constant adiabatic index
GAMMA = 5./3.
# Initialise the Riemann Solver
solver=RiemannSolver(GAMMA)
# The constant time step
timestep = 0.001
# Set up the cells
cells = []

# Loop for 100 cells
# Calculate initial values of the attributes for each cell
for i in range(100):
    cell = Cell() 
    cell._volume = 0.01
    if i < 50:
        #This is the high pressure and density half
        cell._mass = 0.01
        cell._energy = 0.01 / (GAMMA - 1.)
    else:
        #This is the low density half
        cell._mass = 0.00125
        cell._energy = 0.001 / (GAMMA - 1.)
    cell._momentum = 0.
    cell._midpoint = (i+0.5)/100.
    # set the neighbour of the previous cell
    #I am uncertain about this part
    cells.append(cell)
    #cells[-1]._right_ngb = cell

    cell._right_ngb=cells[i-1]

cells[0]._right_ngb=cells[0]

print([cell._right_ngb._mass for cell in cells[49:51]],[cell._mass for cell in cells[49:51]])

i=0.1
while i<0.2:
    for cell in cells:
        # Set up the primative variables
        volume = cell._volume
        mass = cell._mass
        momentum = cell._momentum   
        energy = cell._energy
        density = mass / volume
        velocity = momentum / mass  
        pressure = (GAMMA - 1) * (energy / volume - 0.5 * density * velocity * velocity)    
        
        cell._density = density
        cell._velocity = velocity   
        cell._pressure = pressure   
    
    for cell in cells:
        # Do the boundary conditions
        cell_right = cell._right_ngb        
        densityL = cell._density
        velocityL = cell._velocity
        pressureL = cell._pressure
        densityR = cell_right._density
        velocityR = cell_right._velocity
        pressureR = cell_right._pressure
        densitysol,velocitysol,pressuresol,flag = solver.solve(densityL, velocityL, pressureL, densityR, velocityR, pressureR)
        
        flux_mass = densitysol * velocitysol
        flux_momentum = densitysol * velocitysol * velocitysol + pressuresol
        flux_energy = (((pressuresol*GAMMA) / ((GAMMA - 1))) + (0.5 * densitysol* velocitysol * velocitysol))*velocitysol
    
        
        A = cell._surface_area
        cell._mass = cell._mass - flux_mass * A * timestep
        cell._momentum = cell._momentum - flux_momentum * A * timestep
        cell._energy = cell._energy - flux_energy * A * timestep
        
        cell_right._mass = cell_right._mass + flux_mass * A * timestep
        cell_right._momentum = cell_right._momentum + flux_momentum * A * timestep
        cell_right._energy = cell_right._energy + flux_energy * A * timestep
    
    i = i + timestep
    
plt.plot([ cell._midpoint for cell in cells], [cell._density for cell in cells])
