#!/usr/bin/python
# -*- coding:utf-8 -*-

# Process the experiment files for the R program and creates a corrected file in the current folder
# Uses the bash script for the linebreaks and alias then delete the intermediate file

import sys
import os
import re

if (len(sys.argv) <= 1):
	print("You must indicate at least one file to be processed")
	sys.exit(2)

for id in range(1, len(sys.argv)):
	name = sys.argv[id].split("/")[-1] # Remove the absolute path to place the processed file in the current directory
	os.system("cp " + sys.argv[id] + " ./")
	#name = sys.argv[id]
	os.system("correct_input_files.sh " + name) # Bash script for alias and linebreak

	with open(name + ".corrected") as file:
		nname = name[0:re.search(".txt", name).start()] + ".data" # Changes the extension
		new = open(nname, "w")
		for line in file:
			if (line.split("\t")[2] != "m" and line.split("\t")[2] != "M"):
				new.write(line)
		new.close()
		
# Display the file content to check if everything went fine
		new = open(nname, "r")
		for line in new:
			print(line)
		new.close()

#Deletes the bash output and the copied input
	os.system("rm " + name + ".corrected")
	os.system("rm " + sys.argv[id])


