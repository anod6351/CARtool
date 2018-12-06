def changestart(filepath):

	with open('filepath', 'r') as regions_file:
		with open('regions_fixed_file.bed', 'w') as file:
			for line in regions_file:
				element=line.split('\t')
				row = str(element[0]) + '\t' + str(int(element[1])-1) + '\t' + str(element[2]) + '\t' + str(element[3])
				file.write(row)
	regions_file.close()
	file.close()
	return '/regions_fixed_file.bed'
