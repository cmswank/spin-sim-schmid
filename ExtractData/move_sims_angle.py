import os
import sys
import shutil
import time

fd = "/home/ezra/spin-sim-SD/"
frn = 23000 	# first run number

alpha_steps = 16
beta_steps = 20
gamma_steps = 8

nor = alpha_steps*beta_steps*gamma_steps*2

for rn in range(frn,frn+nor):
	name = "run" + str(rn) + ".root"
	shutil.move(fd + name, fd+ "ExtractData/"+name)

for rn in range(frn,frn+nor):
	cmd = "./extractrootbinary " + str(rn)
	os.system(cmd)
	time.sleep(0.3)


