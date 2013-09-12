#!/Users/work/Documents/Github/work_python/bin/python

### Overview
#This program takes the results of knife edge measurements, which should have the form of 
#an error function, and turns them into a 1/e^2 beam size.
#The short cut is to take only the 90% of max power and the 10% of max power positions.

import numpy as np #import numpy functions
import matplotlib.pyplot as plt #import matplotlib functions
import csv #import csv functions

### Options
fontsize = 16;
shortcut1090 = 1; #%1 if you want to use positions of only 90% and 10% maximum power to determine beam size (shortcut for beam profiling)

if shortcut1090 == 1:
	###Import data
	directoryname = '/Users/work/Documents/Nanofibers/Laser Profiling/Probe Beam Profiling/After Nanofiber/'; 
	inputfilename = 'ProbeLaser_NanofiberOutput_Telescope1_20130906.dat';
	outputfilename = 'ProbeLaser_NanofiberOutput_Telescope1_20130906.txt';
	filename = directoryname + inputfilename;
	
	measurementposition = []; position90=[]; position10=[]; #initialize vectors
	inputfile = open(filename, 'r') #open file for reading
	for line in inputfile:
		line = line.strip() #read in a line at a time
		if not line.startswith("%"): #only read a line if it doesn't start with "%", the comment symbol I've been using in Matlab.
			columns = line.split()
			measurementposition = measurementposition + [float(columns[0])]; #append value of measurementposition from this line to the list of measurementposition values
			position90 = position90 + [float(columns[1])];
			position10 = position10 + [float(columns[2])];
	inputfile.close()

	### Convert lists into numeric vectors
	measurementposition = np.array(measurementposition); #[cm]
	position90 = np.array(position90); #[um]
	position10 = np.array(position10); #[um]

	### Calculate waist size based on Siegman formula (see p. 94 Pascal's Lab Notebook 1)
	waistsize = abs(position90-position10)/1.28; #[um] 1/e^2 beam radius as a function of measurement position along length of beam
    
	### Plot waist size versus distance from lens
	plt.scatter(measurementposition,waistsize,'ro'); #,'sr','MarkerFaceColor','r'
	plt.show()
    #hold on
    #set(gca,'FontSize',fontsize,'FontWeight','bold');
    #xlabel('Longitudinal Position [cm]','FontSize',fontsize,'FontWeight','bold');
    #ylabel('Waist Size [mm]','FontSize',fontsize,'FontWeight','bold');
    #%xlim([0 1.13]);
    #ylim([0 max(waistsize(:))+0.1.*max(waistsize(:))]);
    #%legend('Vertical Size','Horizontal Size','Location','best');
    #%title('Beam Profile of 2988 Laser','FontSize',fontsize,'FontWeight','bold');

	### Manually test
	print waistsize
	print 'status = done'

### Write data to output file
outputdata = np.column_stack((measurementposition,waistsize)); #[cm mm] data is output 

outputfilepath = directoryname + outputfilename;
outputfile = open(outputfilepath,'wb');
filewriter = csv.writer(outputfile, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)

for row in outputdata:
	filewriter.writerow(row)
	
outputfile.close()