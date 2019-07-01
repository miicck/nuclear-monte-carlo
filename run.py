#!/usr/bin/python
import os

# Gather list of directories that need running
dirs = []
for d in os.listdir("samples"):
	if not d.startswith("sample"): continue
	d = "samples/"+d
	if not os.path.isdir(d): continue
	dirs.append(d)

# Run calculations for each sample
for i, d in enumerate(dirs):
	print "Running {0} ({1}/{2})".format(d,i+1,len(dirs))
	os.system("cd "+d+"; nice -15 mpirun pw.x <scf.in> scf.out")
	os.system("cd "+d+"; nice -15 mpirun dos.x <dos.in> dos.out")
