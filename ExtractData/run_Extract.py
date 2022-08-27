import os
import math
import random
import sys
import shutil
import time

#the first run number must be even. (I don't think this is true). 


if len(sys.argv) > 1:
	frn = int(sys.argv[1])	# first run number, must be even
if len(sys.argv) > 2:
	nors = int(sys.argv[2])	# number of run sets

for rn in range(0,nors):
	
	sim_extract = "./extractrootbinary " + str(rn+frn)
	#os.system("gnome-terminal -e 'bash -c \"" +sim_run+ " ;exec bash \"'")
	os.system("nohup " + sim_extract + " > /dev/null 2>&1 &")
	print('extracting run ' + str(rn+frn))
	time.sleep(1)




