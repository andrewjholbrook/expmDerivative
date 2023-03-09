# This file generate the xml files to run BEAST analyses
# Xiang Ji
import argparse, os, re

def main():
	with open("expmDerivative/code/covid10Mar_285plusNYC.Ed.xml",'r') as infile:
		with open("expmDerivative/data/locations.txt", 'w') as outfile:
			for s in infile:
				if re.match("(.*)<state code=\"(.*)", s):
					result = re.search('<state code=\"(.*)\"/>', s)
					outfile.write(str(result.group(1)) + "\n")

if __name__ == '__main__':
	main()

# python get_locations_times.py
