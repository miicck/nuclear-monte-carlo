#!/usr/bin/python
import os
import numpy as np
import matplotlib.pyplot as plt

for d in os.listdir("samples"):

	if not d.startswith("sample"): continue
	d = "samples/"+d
	if not os.path.isdir(d): continue

	# Parse the electron PDOS
	es = None
	pdos = []
	labels = []
	sort_by = []
	for f in os.listdir(d):

		if not f.startswith("pwscf.pdos"): continue
		if not "atm#" in f: continue
		if not "wfc#" in f: continue

		# Get atom number and wfc number from filename
		atm_num = int(f.split("#")[1].split("(")[0])
		wfc_num = int(f.split("#")[2].split("(")[0])

		# Read lines of file in
		f = d+"/"+f
		f = open(f)
		lines = f.read().split("\n")
		f.close()

		# Parse E, ldos(?), pdos from each line
		# (not including first and last line)
		data = []
		for l in lines[1:-1]: data.append([float(w) for w in l.split()])
		data = np.array(data).T
		es   = data[0]
		pdos.append(data[2])

		# Label pdos by atom and wavefunction numbers
		labels.append("Atom {0}, wfc {1}".format(atm_num, wfc_num))

		# Sort pdos by wfc number, then atom number
		# (assumes < 1024 atoms)
		sort_by.append(wfc_num*1024 + atm_num)

	# Sort pdos's according to entries of sort_by
	sort_by, pdos, labels = zip(*sorted(zip(sort_by, pdos, labels)))

	# Plot pdos
	total = np.zeros(len(pdos[0]))
	for p, l in zip(pdos, labels):
		plt.fill_between(es, total, total+p, label=l)
		total += p
	plt.legend()
	plt.show()
