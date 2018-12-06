def detail_samtools(Regions, Read_depth):

	# create a detailed list with all depth values from the same region in a sub list. From samtools depth calculations
	# samtools generates a depth file with: chr, position and coverage depth value 
	# Regeions comes from the bed file with chr, start, stop, region name

	detailed =[]
	list_temp=[]
	previous_chr = Read_depth[0][0]
	Region_row = 0
	count=0
	index = 0
	size_list = len(Read_depth)

	for line in Read_depth:
		Region_row = Regions[index]
		if str(line[0]) == str(previous_chr) and (int(line[1])) <= int(Region_row[2]):
			list_temp.append(line[2])

		else:
			previous_chr=line[0]
			detailed.append(list_temp)
			list_temp=[]
			list_temp.append(line[2])
			index+=1

		count+=1
		if count == size_list:
			detailed.append(list_temp)

	return detailed



