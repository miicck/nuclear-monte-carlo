#!/usr/bin/python
import numpy as np
import sys
import os

# Check if we want to plot the sampled nuclear density
plot = "plot" in sys.argv[1:]
if plot:
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib.pyplot as plt

def nuclear_density(x):
	return np.exp(-np.linalg.norm(x)**2)

# Sample the function metro according to the metropolis algorithm
#       x0       = start point
# 	steps    = # of metropolis iterations
#       samples  = # of samples to take from those metropolis steps
#	metro    = function to sample
# returns a pair [weights, path] where weights are the weights
# corresponding to each sampled point in path.
def generate_metro_path(x0, 
			steps, 
                        samples,
			metro=nuclear_density):

	# Start path at initial position x0 with weight 1.0
	path    = [x0]
	weights = [1.0]

	# Sample steps using metropolis algorithm
	# scales step size dynamically
	step_size = 0.1
	accepted  = 0.0
	rejected  = 0.0

	# We do steps-1 iterations so when we include x0
	# the total weight = steps
	for step in range(0, steps-1):

		# Make a random step and clamp to fractional coordinates
		xnew = path[-1] + np.random.normal(0, step_size, path[-1].shape)
		xnew = np.mod(xnew, 1)

		# Evaluate the metropolis function
		# before and after the step
		m1 = metro(path[-1])
		m2 = metro(xnew)

		# Accept step with probability m2/m1
		if m2/m1 > np.random.rand():
			# Add new point with weight = 1.0
			path.append(xnew)
			weights.append(1.0)
			accepted += 1.0
		else:
			# Stay where I am (increase weight by 1.0)
			weights[-1] += 1.0
			rejected += 1.0

		# Modify step size to get close to optimal acceptance ratio
		# of 0.234 (see The Annals of Applied Probability
		# Vol. 7, No. 1 (Feb., 1997), pp. 110-120)
		acceptance_ratio = accepted/(accepted+rejected)
		step_size += 0.1*(acceptance_ratio-0.234)
		if step_size < 0: step_size = 0

        # Resample to avoid serial correlation
        new_weights = []
        new_path    = []
        for j in range(0, len(weights), int(len(weights)/samples)):
            if len(new_weights) >= samples: break
            new_weights.append(weights[j])
            new_path.append(path[j])
        
	print "Acceptance ratio: ", acceptance_ratio*100, "% (optimal = 23.4%)"
        print "Requested samples: {0}, produced samples: {1}".format(samples, len(new_weights))
	return [np.array(new_weights), np.array(new_path)]

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

	# Convert to np array for convinience
	positions = np.array(positions)

	# Check positions are fractional
	for c in positions.flatten():
		if c >= 0 and c <= 1.0: continue
		print "Error, expected fractional coordinates in [0,1], but got a coordinate: ", c

	return positions

# Creates the directory and input files for a given weight/position
direc_count = 0
def create_input(weight, positions, to_copy="scf.in"):
	global direc_count

	# Create the directory where the input will be put
	direc_count += 1
	direc = "samples/sample_"+str(direc_count)
	os.system("mkdir samples 2>/dev/null")
	os.system("mkdir "+direc)
	os.system("cp dos.in "+direc)
	os.system("echo "+str(weight)+" > "+direc+"/weight")

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

# A guassian centred around each initial position
def guassian_spread(positions, init_positions, decay_length=0.1):

	# Work out the coordinates of the guassian centres
	# up to some range
	centres = []
	max_range = 1.0
	lattice = np.arange(-max_range, max_range+0.5, 1.0)
	for i in init_positions:
		for dx in lattice:
			for dy in lattice:
				for dz in lattice:
					centres.append(i + np.array([dx,dy,dz]))

	result = 1.0
	for p in positions:

		# Sum up the probability of p to be where it
		# is given the centres
		vp = 0
		for c in centres:
			r = np.linalg.norm(p-c)/decay_length
			vp += np.exp(-r**2)

		# Prob(r1,r2,r3..) = prob(r1)*prob(r2)*prob(r3)...
		result *= vp
			
	return result

# Don't overwrite previous run
if os.path.isdir("samples"):
	if "-f" in sys.argv:
		os.system("rm -r samples")
	else:
		print "Error, refusing to overwirte samples directory."
		quit()

# Check we've got the right number of arguments
if len(sys.argv) < 2:
	print "Error, I require one argument, the number of MC samples!"
	quit()

# Get number of samples
try:
	samples = int(sys.argv[1])
except:
	print "Could not parse number of samples from ", sys.argv[1]
	quit()

init_pos = read_positions()
guass_dist = lambda x : guassian_spread(x, init_pos)
weights, path = generate_metro_path(init_pos, samples*10, samples, metro=guass_dist)
for w, p in zip(weights, path): create_input(w, p)

if plot:
	# Plot the nuclear density
	ax = plt.gcf().add_subplot(111, projection='3d')
	ax.scatter(init_pos.T[0], init_pos.T[1], init_pos.T[2])
	ax.scatter(path.T[0], path.T[1], path.T[2])
	plt.show()
