import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-f","--filepath", help="Write the path to the Bed file, changes start position from 1-based to 0-based")
args = parser.parse_args()

with open(args.filepath, 'r') as regions_file:
	with open('regions_fixed_file.bed', 'w') as file:
		for line in regions_file:
			element=line.split('\t')
			row = str(element[0]) + '\t' + str(int(element[1])-1) + '\t' + str(element[2]) + '\t' + str(element[3])
			file.write(row)