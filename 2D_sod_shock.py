# Fluids 1D Sod shock

from riemannsolver import RiemannSolver
import matplotlib.pyplot as plt
import sys
import numpy as np
import copy as copy

###############################################################################

# The Cell class
# This sets up the parameters of the cell
class Cell(object):
    def __init__(self):
        self._volume = 0.
        self._mass = 0.
        self._xmomentum = 0.
        self._ymomentum = 0.
        self._energy = 0.
        self._density = 0.
        self._xvelocity = 0.
        self._yvelocity = 0.
        self._xpressure = 0.
        self._ypressure = 0.
        self._right_ngb = None
        self._up_ngb = None
        self._right_area = 1.
        self._up_area = 1.
 
###############################################################################
       
# The constant adiabatic index
GAMMA = 5./3.
# Initialise the Riemann Solver
solver = RiemannSolver(GAMMA)
# Set up the cells
ncell = 100
cells = [ncell,ncell]
xnormal = 1.
ynormal = 1.
flag=-1

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
    cell._xmomentum = 0.*xnormal
    cell._ymomentum = 0.*ynormal
    cell._xmidpoint = (i+0.5)/ncell
    cell._ymidpoint = (i+0.5)/ncell
    # Add cell to list cells
    cells[i:].append(cell)
    cells[:i].append(cell)
    # Set up right neightbour of cells
    cell._right_ngb=cells[(i-1):]
    cell._up_ngb=cells[:(i-1)]

###############################################################################

def xflux(cell):

    velocityL = cell._xvelocity*xnormal
    densityL = cell._density
    pressureL = cell._xpressure        
    cell_right = cell._right_ngb
    velocityR = cell_right._xvelocity*xnormal
    densityR = cell_right._density
    pressureR = cell_right._xpressure
    #Riemann solver
    densitysol,velocitysol,pressuresol,flag = solver.solve(densityL, velocityL, pressureL, densityR, velocityR, pressureR)

    #Calculate the fluxes
    flux_mass = densitysol * velocitysol
    flux_momentum = ((densitysol * velocitysol * velocitysol) + pressuresol) *xnormal
    flux_energy = (((pressuresol*GAMMA) / (GAMMA - 1)) + (0.5 * densitysol* velocitysol * velocitysol))*velocitysol

    A = cell._right_area
    
    #Update cell parameters
    cell._mass = cell._mass - flux_mass * A * timestep
    cell._xmomentum = (cell._xmomentum - flux_momentum * A * timestep)*xnormal
    cell._energy = cell._energy - flux_energy * A * timestep
    
    #Update right neighbour parameters
    cell_right._mass = cell_right._mass + flux_mass * A * timestep
    cell_right._xmomentum = (cell_right._xmomentum + flux_momentum * A * timestep)*xnormal
    cell_right._energy = cell_right._energy + flux_energy * A * timestep
    
###############################################################################
def yflux(cell,flag):

    velocityL = cell._yvelocity*ynormal
    densityL = cell._density
    pressureL = cell._ypressure        
    cell_right = cell._up_ngb
    velocityR = cell_right._yvelocity*ynormal
    densityR = cell_right._density
    pressureR = cell_right._ypressure
    
    if (flag==-1):
        velocityL=velocityL
    elif flag==1:
        velocityL=velocityR
    
    #Riemann solver
    densitysol,velocitysol,pressuresol,flag = solver.solve(densityL, velocityL, pressureL, densityR, velocityR, pressureR)

    #Calculate the fluxes
    flux_mass = densitysol * velocitysol
    flux_momentum = ((densitysol * velocitysol * velocitysol) + pressuresol)*ynormal
    flux_energy = (((pressuresol*GAMMA) / (GAMMA - 1)) + (0.5 * densitysol* velocitysol * velocitysol))*velocitysol

    A = cell._up_area
    
    #Update cell parameters
    cell._mass = cell._mass - flux_mass * A * timestep
    cell._ymomentum = (cell._ymomentum - flux_momentum * A * timestep)*ynormal
    cell._energy = cell._energy - flux_energy * A * timestep
    
    #Update right neighbour parameters
    cell_right._mass = cell_right._mass + flux_mass * A * timestep
    cell_right._ymomentum = (cell_right._ymomentum + flux_momentum * A * timestep)*ynormal
    cell_right._energy = cell_right._energy + flux_energy * A * timestep
    
###############################################################################

i=0.0               #Initial time
timestep = 0.001    #The constant time step
total_time=0.07

# Loop over time
while i<total_time:

    for cell in cells:
        # Set up the primative variables
        volume = cell._volume
        mass = cell._mass
        xmomentum = cell._xmomentum * xnormal  
        ymomentum = cell._ymomentum *ynormal
        energy = cell._energy
        
        if (volume==0 or mass==0):
            print("ERROR: volume or mass is equal to zero \n",
                  "Mass = ", mass, "\n Volume = ", volume,)
            sys.exit()
        
        density = mass / volume
        xvelocity = (xmomentum / mass)*xnormal
        yvelocity = (ymomentum / mass)*ynormal
        xpressure = (GAMMA - 1) * (energy / volume - 0.5 * density * xvelocity * xvelocity)    
        ypressure = (GAMMA - 1) * (energy / volume - 0.5 * density * yvelocity * yvelocity)    
        
        if(mass<=0 or density<=0 or energy<=0 or xpressure<=0 or ypressure<=0):
            print("ERROR: mass, density, energy or pressure is less than or equal to zero \n", 
                  "Mass = ", mass, "\n Density = ", density, "\n Energy = ", energy, "\n xPressure = ", xpressure, "\n yPressure = ", ypressure)
            sys.exit()

        cell._density = density
        cell._xvelocity = xvelocity*xnormal 
        cell._yvelocity = yvelocity*ynormal 
        cell._xpressure = xpressure   
        cell._ypressure = ypressure   
        
    cells[0]._right_ngb=cells[99]
    cells[0]._up_ngb=cells[99]
    
    for cell in cells:
        
        xflux(cell)
        yflux(cell,flag)

    i = i + timestep    #Increase the time 

###############################################################################
plt.scatter([cell._xmidpoint for cell in cells], [cell._density for cell in cells], label="Density",s=5)
plt.scatter([cell._ymidpoint for cell in cells], [cell._density for cell in cells], label="Density",s=5)
#plt.scatter([cell._xmidpoint for cell in cells], [cell._xvelocity for cell in cells], label="Density",s=5)
#plt.scatter([cell._ymidpoint for cell in cells], [cell._yvelocity for cell in cells], label="Density",s=5)
plt.legend()
plt.xlabel("Cell Position")
plt.ylabel("Cell Density") 
