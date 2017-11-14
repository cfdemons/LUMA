# Script that changes the restart file from LUMA 1.6.1 or older to the format that allows restarting with a differnt value of dt (LUMA versions after LUMA 1.6.1). 
# It only works with a single grid, not multi-grid. 
# It only works with D3Q19. 
# Author: Marta Camps Santasmasas
# Date: 08/11/2017

import numpy as np

# Input folder
finput = './/restartFiles203flow05dt35cellsdP00589//'

# Output folder
foutput = './/input//'

# Number of dimensions
dims = 3

# Number of processes
nproc = 64

# Discretization
dx = 1/35.0
dt = 0.0016

# Relaxation frequency
omega = 1.0/0.501535

# Data for feq
c_opt = np.matrix( [ [ 1, 0, 0 ], [ -1, 0, 0 ], [ 0, 1, 0 ], [ 0, -1, 0 ],[ 0, 0, 1 ],[ 0, 0, -1 ],[1, 1, 0 ],[ -1, -1, 0 ],[ 1, -1, 0 ],[ -1, 1, 0 ],[ 0, 1, 1 ],[ 0, -1, -1 ],[ 0, 1, -1 ],[ 0, -1, 1 ],[1, 0, 1 ],[ -1, 0, -1 ],[ -1, 0, 1 ],[ 1, 0, -1 ],[0, 0, 0, ]])
cs = 1.0/np.sqrt(3.0)
w = np.array([1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0,1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,1.0/3.0])

# Writing format
format = '%1d\t%1d'
for i in range(0,(23+dims)):
	format = format + '\t%1.8e'

def SQ(num):
	return num*num
	
def f_equilibrium(ux, uy, uz, rho, v, dims):
	# Compute the parts of the expansion for feq
	
	if dims == 3:
		A = (c_opt[v,0] * ux) +  (c_opt[v,1] * uy) + (c_opt[v,2] * uz)
		B = (SQ(c_opt[v,0]) - SQ(cs)) * SQ(ux) + (SQ(c_opt[v,1]) - SQ(cs)) * SQ(uy) + (SQ(c_opt[v,2]) - SQ(cs)) * SQ(uz) +  2 * c_opt[v,0] * c_opt[v,1] * ux * uy + 2 * c_opt[v,0] * c_opt[v,2] * ux * uz + 2 * c_opt[v,1] * c_opt[v,2] * uy * uz
	elif dims == 2:
		A = (c_opt[v,0] * ux) + (c_opt[v,1] * uy)
		B = (SQ(c_opt[v,0]) - SQ(cs)) * SQ(ux) + (SQ(c_opt[v,1]) - SQ(cs)) * SQ(uy) + 2 * c_opt[v,0] * c_opt[v,1] * ux * uy
	else:
		print('Wrong number of dimensions')
		return 0
		
	# Compute f^eq
	fequi = rho * w[v] * ( 1.0 + (A / SQ(cs)) + (B / (2.0 * SQ(cs)*SQ(cs)) ) ) 
	
	return fequi
	


# Loop over all files 
for f in range(0,nproc):
	
	# Load data. 
	name ='restart_LBM_Rnk' + str(f) + '.out'
	print(name)
	oldRestart = np.loadtxt(finput + name, delimiter='\t',unpack=True)
	
	# Extract data
	ux = oldRestart[24,:]
	uy = oldRestart[25,:]
	uz = oldRestart[26,:]
	rho = oldRestart[24+dims,:]
	
	# Create new matrix to store data by copying the old matrix. 
	newRestart = np.copy(oldRestart)
	
	# Copy data in the correct order and modify it as needed
	# density
	newRestart[(5+dims),:] = rho
	
	# velocity
	newRestart[5,:] = (ux*dx)/dt
	newRestart[6,:] = (uy*dx)/dt
 	newRestart[7,:] = (uz*dx)/dt
		
	# Non-equilibrium f
	for v in range(0,19):
		f = oldRestart[(5+v),:]
		feq = f_equilibrium(ux,uy,uz,rho,v, dims)
		fneq = ((f-feq)*omega)/(feq*dt)
		newRestart[6+dims+v,:] = fneq

	# write data to new file
	np.savetxt(foutput + name, newRestart.transpose(), delimiter='\t', fmt=format)
	













