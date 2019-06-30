from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

def nuclear_density(x):
	return np.exp(-np.linalg.norm(x)**2)

# Sample the function metro according to the metropolis algorithm
#       x0       = start point
# 	steps    = # of metropolis iterations
# 	max_step = maximum step size in config space
#	metro    = function to sample
# returns a pair [weights, path] where weights are the weights
# corresponding to each sampled point in path
def generate_metro_path(x0, 
			steps=1000, 
			max_step=0.01, 
			metro=nuclear_density):

	# Start path at initial position x0 with weight 1.0
	path    = [x0]
	weights = [1.0]

	# Sample steps using metropolis algorithm
	for step in range(0, steps):

		# Make a random step
		delta = 1-2*np.random.rand(*path[-1].shape)
		delta *= max_step
		xnew = path[-1] + delta

		# Evaluate the metropolis function
		# before and after the step
		m1 = metro(path[-1])
		m2 = metro(xnew)

		# Accept step with probability m2/m1
		if m2/m1 > np.random.rand():
			# Add new point with weight = 1.0
			path.append(xnew)
			weights.append(1.0)
		else:
			# Stay where I am (increase weight by 1.0)
			weights[-1] += 1.0

	return [np.array(weights), np.array(path)]

# Read the positions of the atoms in the q.e ".in" file
def read_positions(filename="scf.in"):
	
	# Read lines
	f = open(filename)
	lines = f.read().split("\n")
	f.close()

	# Parse positions
	positions = None
	for l in lines:
		l = l.lower()
		if "atomic_positions" in l:
			# We've hit the atomic_positions
			# label, we can start parsing them now
			positions = []
			continue

		if positions is None: 
			# Havent hit atomic_positions
			# label yet
			continue

		try:
			# Attempt to parse atom name, x, y, and z
			n,x,y,z = l.split()
			positions.append([float(c) for c in [x,y,z]])
		except:
			# Parsing failed => we've run out of atoms
			break

	return np.array(positions)

# Creates the directory and input files for a given weight/position
direc_count = 0
def create_input(weight, positions, to_copy="scf.in"):
	global direc_count

	# Create the directory where the input will be put
	direc_count += 1
	direc = "samples/sample_"+str(direc_count)
	os.system("mkdir samples 2>/dev/null")
	os.system("mkdir "+direc)

	# Copy the scf input file there, modifying it
	# to contin the correct atom positions
	fw = open(direc+"/"+to_copy, "w")
	fr = open(to_copy)
	lines = fr.read().split("\n")
	fr.close()

	atom_index = -1
	for l in lines:	

		# Replace atom coordinates
		if atom_index >= 0 and atom_index < len(positions):
			l = l.split()[0]+"".join([" "+str(c) for c in positions[atom_index]])
			atom_index += 1

		# Flag to start replacing atom coordinates
		if "atomic_positions" in l.lower():
			atom_index = 0

		# Copy the (potentially modified) line
		fw.write(l+"\n")

	fw.close()

def guassian_spread(positions, init_positions):
	exps = [np.exp(-np.linalg.norm(p-i)**2) for p,i in zip(positions, init_positions)]
	return np.product(exps)

# Don't overwrite previous run
if os.path.isdir("samples"):
	print "Error, refusing to overwirte samples directory."
	quit()


init_pos = read_positions()
guass_dist = lambda x : guassian_spread(x, init_pos)

weights, path = generate_metro_path(init_pos, metro=guass_dist)
create_input(weights[0], path[-1])


ax = plt.gcf().add_subplot(111, projection='3d')
ax.scatter(path.T[0], path.T[1], path.T[2], alpha=0.1)
plt.show()
