# Fluids 1D Sod shock

from riemannsolver import RiemannSolver
import matplotlib.pyplot as plt
import copy

###############################################################################

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
 
###############################################################################
       
# The constant adiabatic index
GAMMA = 5./3.
# Initialise the Riemann Solver
solver=RiemannSolver(GAMMA)
# Set up the cells
cells = []

###############################################################################

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
    # Add cell to list cells
    cells.append(cell)
    # Set up right neightbour of cells
    cell._right_ngb=cells[i-1]

###############################################################################
'''
#Setting up initial consitions and plotting
for cell in cells:
    # Set up the primative variables
    volume = cell._volume
    mass = cell._mass
    momentum = cell._momentum   
    energy = cell._energy
    density = mass / volume
    velocity = momentum / mass  
    pressure = (GAMMA - 1) * ((energy / volume) - (0.5 * density * velocity * velocity)) 
    print(velocity)
    
    cell._density = density
    cell._velocity = velocity   
    cell._pressure = pressure 

plt.plot([cell._midpoint for cell in cells], [cell._density for cell in cells], label="Density")
plt.plot([cell._midpoint for cell in cells], [cell._velocity for cell in cells], label="Velocity")
plt.plot([cell._midpoint for cell in cells], [cell._pressure for cell in cells], label ="Pressure")
plt.title("Initial Conditions")
plt.xlabel("Cell Position")
plt.legend()
'''
###############################################################################
'''
#Storing the first and last cell values
left_bound=cells[0]
left_bound._velocity=-cells[0]._velocity
right_bound=cells[(len(cells)-1)]
right_bound._velocity=-cells[(len(cells)-1)]._velocity

# Checking that the properties change in the middle for initial consitions
#print([cell._right_ngb._mass for cell in cells[49:51]],[cell._mass for cell in cells[49:51]])
'''
###############################################################################

def flux(cell, velocityL, velocityR, densityR, pressureR):

    cell_right = cell._right_ngb 

    densityL = cell._density
    pressureL = cell._pressure
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

j=0                 #Counter
i=0.0               #Initial time
k=0
timestep = 0.001    #The constant time step


# Loop over time
while i<0.2:
    k=0

    for cell in cells:
        # Set up the primative variables
        
        volume = cell._volume
        mass = cell._mass
        #print(j,k,mass)
        momentum = cell._momentum   
        energy = cell._energy
        density = mass / volume
        velocity = momentum / mass  
        pressure = (GAMMA - 1) * (energy / volume - 0.5 * density * velocity * velocity)    
        
        cell._density = density
        cell._velocity = velocity   
        cell._pressure = pressure   
        k=k+1
        
    Lcell=copy.deepcopy(cells[0])
    Lcell._right_ngb=cells[0]

    
    Rcell=copy.copy(cells[99])
    #Rcell._velocity=-Rcell._velocity
    cells[99]._right_ngb=Rcell

    
    #flux(cell, velocityL, velocityR, densityR, pressureR)

    
    for cell in cells:
        
        # Deal with the ghost neighbours at the boundaries
        if cell==cells[0]:
            #The leftmost cell
            cell_right = Lcell
            velocityL = -cell._velocity
            velocityR = cell_right._velocity
            densityR = cell_right._density
            pressureR = cell_right._pressure
            #print(j,"left",velocityL,velocityR)
            
        elif cell==cells[99]:
            #The rightmost cell
            cell_right = Rcell            
            velocityL = cell._velocity
            velocityR = -cell_right._velocity 
            densityR = cell_right._density
            pressureR = cell_right._pressure
            #print(j," right ",velocityL,velocityR)
            
        else:
        
            #The middle cells
            cell_right = cell._right_ngb 
            velocityL = cell._velocity
            
            velocityR = cell_right._velocity
            densityR = cell_right._density
            pressureR = cell_right._pressure
            #print(j," middle ", velocityL, velocityR)
        
        #Left cell values
        densityL = cell._density
        pressureL = cell._pressure

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
        
        #if (j>=90):
            #print(j,cell_right._mass)
        
        j=j+1           #Counter
    i = i + timestep    #Increase the time 
    #plt.plot([cell._midpoint for cell in cells], [cell._density for cell in cells], label="Density")

###############################################################################


###############################################################################

#Plot the midpoint of each cell versus the density of the cell    
plt.plot([cell._midpoint for cell in cells], [cell._density for cell in cells], label="Density")
#plt.plot([cell._midpoint for cell in cells], [cell._velocity for cell in cells], label="Velocity")
#plt.plot([cell._midpoint for cell in cells], [cell._pressure for cell in cells], label ="Pressure")
#plt.legend()

