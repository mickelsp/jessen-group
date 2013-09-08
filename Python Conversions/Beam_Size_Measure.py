#This program takes the results of knife edge measurements, which should have the form of 
#an error function and turns them into a 1/e^2 beam size.

fontsize = 16;
shortcut1090 = 1; #%1 if you want to use positions of only 90% and 10% maximum power to determine beam size (shortcut for beam profiling)

if shortcut1090 == 1:
	#Import data
	filepath = '/Users/work/Documents/Nanofibers/Laser Profiling/Probe Beam Profiling/After Nanofiber/';
	inputfilename = 'ProbeLaser_NanofiberOutput_Telescope1_20130906.dat';
	outputfilename = 'ProbeLaser_NanofiberOutput_Telescope1_20130906.txt';
	filename = filepath + inputfilename;
	
	f = open(filename, 'r')
	for line in f:
		line = line.strip()
		if not line.startswith("%"):
			columns = line.split()
	f.close()

#Manually test
#print line
print 'status = done'
