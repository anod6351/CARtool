def detail(Regions, Reads):

	import sys
	import os
	####################### Computations of per base/amplicon coverage using bedtools ##############################
	# Compute the per base coverage for the regions in the bedfile, columns in Region_coverage: chromosome, region start, region stop, name, index number of base (what base in the region), coverage
	
	command = "bedtools coverage -a " + str(Regions) + " -b " + str(Reads) + " -d > CAR_output/detailedCoverage.bed"
	os.popen(command)

	# Create per position coverage for the strand specific bam files 
	detailed_file = open("CAR_output/detailedCoverage.bed", "r")
	Region_coverage_perbase= []

	for line in detailed_file:
		Region_coverage_perbase.append(line)
	#print(Region_coverage_perbase[0])



	# create a detailed list
#Region_coverage_perbase = [['ch1',0,3, 'Gene1', 1, 20], ['ch1',0,3, 'Gene1', 2, 21], ['ch1',0,3, 'Gene1', 1, 22],['ch2', 0,5, 'Gene2', 1,30], ['ch1',0,3, 'Gene1', 1, 22]]
	detailed =[]
	list_temp=[]
	previous_gene = Region_coverage_perbase[0]
	previous_gene = previous_gene[3]
	print(previous_gene)
	count=0
	size_list = len(Region_coverage_perbase)
	print(size_list)

	for line in Region_coverage_perbase:
		count+=1
		#print(line)

		if str(line[3]) == str(previous_gene):
			list_temp.append(line[5])

		else:
			previous_gene=line[3]
			detailed.append(list_temp)
			list_temp=[]
			list_temp.append(line[5])
			
		if count == size_list:
			detailed.append(list_temp)

	#print('Region bedtools')
	#print(Region_coverage_perbase)
	#print(Region_coverage_perbase[0])
	print(Region_coverage_perbase[89])
	print(Region_coverage_perbase[90])


	print('Region bedtools lista')
	print(detailed[0])
	#print(detailed[1])
	#print(detailed[2])
	#print(detailed[3])
	#print(detailed[4])
	#print(detailed[1])
	return detailed
