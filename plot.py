#!/usr/bin/python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Return a weighted standard deviation
def std(vals, ws):
	av  = np.average(vals, weights=ws, axis=0)
	var = np.average((vals-av)**2, weights=ws, axis=0)
	return np.sqrt(var)

directories = os.listdir("samples")
directories = ["samples/"+d for d in directories]
if len(sys.argv) > 1: directories = sys.argv[1:]

weights  = []
energies = []
doss     = []
max_min  = -np.inf
min_max  = np.inf
for d in directories:

	if not os.path.isdir(d):
		continue

	f = d+"/weight"
	if not os.path.isfile(f):
		continue

	# Parse the weight
	f = open(f)
	weight = float(f.read())
	f.close()

	# Parse the electron DOS
	# Read lines of file in
	f = d+"/pwscf.dos"
	if not os.path.isfile(f):
		continue

	f = open(f)
	lines = f.read().split("\n")
	f.close()

	# Parse E, dos, Int pdos from each line
	# (not including first and last line)
	data = []
	efermi = float(lines[0].split("=")[1].split("e")[0])
	for l in lines[1:-1]: data.append([float(w) for w in l.split()])
	data = np.array(data).T
	es   = data[0] - efermi
	dos  = data[1]

	# Find the range of es present in every DOS
	if max(es) < min_max: min_max = max(es)
	if min(es) > max_min: max_min = min(es)
	
	# Record the results
	weights.append(weight)
	energies.append(es)
	doss.append(dos)

# Align all the dos arrays
count = -1
for n_sample in range(0, len(weights)):

	e = energies[n_sample]
	dos = doss[n_sample]

	# Get the indicies bounding our plotted range
	for i in range(0, len(e)):
		if e[i] < max_min:
			i_min = i
		if e[i] > min_max:
			i_max = i
			break

	# Ensure the number of elements in each array is the same
	if count < 0: count  = i_max - i_min
	else: i_max += count - (i_max - i_min)

	# Throw away data outside plotted range
	energies[n_sample] = e[i_min:i_max]
	doss[n_sample]     = dos[i_min:i_max]

# Plot each DOS
plt.subplot(211)
for n_sample in range(0, len(weights)):
	plt.plot(energies[n_sample], doss[n_sample])

sum_weight = np.sum(weights)
sum_dos = np.dot(weights, doss)/sum_weight
sum_energies = np.dot(weights, energies)/sum_weight
std_dos = std(doss, weights)

# Plot the average dos and std deviation
plt.subplot(212)
plt.fill_between(sum_energies, sum_dos-std_dos, sum_dos+std_dos, alpha=0.5, color="black")
plt.plot(sum_energies, sum_dos, color="black", linestyle=":")
plt.show()
