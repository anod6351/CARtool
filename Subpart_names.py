#subpart names

def bar_names_generator(list_to_compress, regions):

#example input, regions = ['hej.g.1.e', 'hej.g.2.e', 'hej.g.1.e', 'h.g.1.e']

	regions_temp_bar = []
	region_name=[]
	region_name_bar=[]


	# compute region subpart name list
	for line in regions:
		l = line.split('.')
		# here the third element in the region name after a '.' is taken as the subpart name 
		region_name.append(l[2])

	# Add subpart names from the same combined region in a nested list
	previous_region = []
	index=-1
	regions.append('end.extra')
	current_region = []
	formated_bar_temp=[]
	detailed_list_formated_bar=[]

						
	for line in regions:
		current_region = line.split('.')
		current_region = current_region[0]

		if not(previous_region ==[]):

			if previous_region == current_region:
				regions_temp_bar.append(region_name[index])
				formated_bar_temp.append(list_to_compress[index])

			else:
				regions_temp_bar.append(region_name[index])
				region_name_bar.append(regions_temp_bar)
				formated_bar_temp.append(list_to_compress[index])
				detailed_list_formated_bar.append(formated_bar_temp)
				
				regions_temp_bar=[]
				formated_bar_temp=[]
					
		index+=1
		previous_region=current_region
	regions.pop()
	
	return region_name_bar, detailed_list_formated_bar



