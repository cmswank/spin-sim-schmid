import shutil

fd = "/home/ezra/spin-sim-SD/ExtractData/"
frn = 23000	# first run number
nor = 5120

for rn in range(frn,frn+nor):
	name = "OutputHeatflush" + str(rn) + ".dat"
	shutil.move(fd + name, fd+ "ang search 4/"+name)