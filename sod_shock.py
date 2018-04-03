# Fluids 1D Sod shock

from riemannsolver import RiemannSolver
import matplotlib.pyplot as plt
import sys
import numpy as np

###############################################################################

# The Cell class
# This sets up the parameters of the cell
class Cell(object):
    def __init__(self):
        self._volume = 0.
        self._mass = 0.
        self._momentum = [0.,0.]
        self._energy = 0.
        self._density = 0.
        self._velocity = [0.,0.]
        self._midpoint = np.zeros((1,1))
        self._pressure = 0.
        self._right_ngb = None
        self._surface_area = 1.
 
###############################################################################
       
# The constant adiabatic index
GAMMA = 5./3.
# Initialise the Riemann Solver
solver=RiemannSolver(GAMMA)
# Set up the cells
cells = []
ncell=100

###############################################################################

# Loop for 100 cells
# Calculate initial values of the attributes for each cell
for i in range(ncell):
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
    cell._midpoint = (i+0.5)/ncell
    # Add cell to list cells
    cells.append(cell)
    # Set up right neightbour of cells
    cell._right_ngb=cells[i-1]

###############################################################################

def flux(cell):

    velocityL = cell._velocity
    densityL = cell._density
    pressureL = cell._pressure        
    cell_right = cell._right_ngb
    velocityR = cell_right._velocity
    densityR = cell_right._density
    pressureR = cell_right._pressure
    #Riemann solver
    densitysol,velocitysol,pressuresol,flag = solver.solve(densityL, velocityL, pressureL, densityR, velocityR, pressureR)

    #Calculate the fluxes
    flux_mass = densitysol * velocitysol
    flux_momentum = (densitysol * velocitysol * velocitysol) + pressuresol
    flux_energy = (((pressuresol*GAMMA) / (GAMMA - 1)) + (0.5 * densitysol* velocitysol * velocitysol))*velocitysol

    A = cell._surface_area
    
    #Update cell parameters
    cell._mass = cell._mass - flux_mass * A * timestep
    cell._momentum = cell._momentum - flux_momentum * A * timestep
    cell._energy = cell._energy - flux_energy * A * timestep
    
    #Update right neighbour parameters
    cell_right._mass = cell_right._mass + flux_mass * A * timestep
    cell_right._momentum = cell_right._momentum + flux_momentum * A * timestep
    cell_right._energy = cell_right._energy + flux_energy * A * timestep
    
###############################################################################

i=0.0               #Initial time
timestep = 0.001    #The constant time step
total_time=0.2
# Loop over time
while i<total_time:

    for cell in cells:
        # Set up the primative variables
        volume = cell._volume
        mass = cell._mass
        momentum = cell._momentum   
        energy = cell._energy
        
        if (volume==0 or mass==0):
            print("ERROR: volume or mass is equal to zero \n",
                  "Mass = ", mass, "\n Volume = ", volume,)
            sys.exit()
        
        density = mass / volume
        velocity = momentum / mass  
        pressure = (GAMMA - 1) * (energy / volume - 0.5 * density * velocity * velocity)    
        
        if(mass<=0 or density<=0 or energy<=0 or pressure<=0):
            print("ERROR: mass, density, energy or pressure is less than or equal to zero \n", 
                  "Mass = ", mass, "\n Density = ", density, "\n Energy = ", energy, "\n Pressure = ", pressure)
            sys.exit()

        cell._density = density
        cell._velocity = velocity   
        cell._pressure = pressure   
        
    cells[0]._right_ngb=cells[99]

    for cell in cells:
        
        flux(cell)

    i = i + timestep    #Increase the time 

###############################################################################
plt.scatter([cell._midpoint for cell in cells], [cell._density for cell in cells], label="Density",s=5)
