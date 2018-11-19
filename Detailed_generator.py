def detail(Regions, Reads):

	import pybedtools

	####################### Computations of per base/amplicon coverage using bedtools ##############################
	# Compute the per base coverage for the regions in the bedfile, columns in Region_coverage: chromosome, region start, region stop, name, index number of base (what base in the region), coverage
	
	Region_coverage_perbase = Regions.coverage(Reads, d=True)

	# create a detailed list

	detailed =[]
	list_temp=[]
	previous_gene = Region_coverage_perbase[0][3]
	print(previous_gene)
	count=0
	size_list = len(Region_coverage_perbase)
	print(size_list)

	for line in Region_coverage_perbase:
		count+=1

		if str(line[3]) == str(previous_gene):
			list_temp.append(line[5])

		else:
			previous_gene=line[3]
			detailed.append(list_temp)
			list_temp=[]
			list_temp.append(line[5])
			
		if count == size_list:
			detailed.append(list_temp)

	return detailed






