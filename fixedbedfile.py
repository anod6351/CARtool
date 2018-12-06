# To change the end position in a bedfile from 0-based to 1-based if needed. CARtool only takes 0-based start positions and 1 based end positions
def changestart(filepath):

	with open(filepath, 'r') as regions_file:
		with open('regions_fixed_file.bed', 'w') as myfile:
			for line in regions_file:
				element=line.split('\t')
				row = str(element[0]) + '\t' + str(element[1]) + '\t' + str(int(element[2])+11) + '\t' + str(element[3])
				myfile.write(row)
	regions_file.close()
	myfile.close()
	return '/regions_fixed_file.bed'
