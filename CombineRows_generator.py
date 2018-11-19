
def CombineRowsList(listToCompress, regions):
	from collections import OrderedDict
	#listToCompress=[[1,2,2,3],[1,2,2,23,4],[2,4,6,7,8], [1,2,3,4]]
	#regions = ['GeneA.g.1.e', 'GeneA.g.2.e', 'GeneA.g.1.e', 'GeneB.g.1.e']

	# Variables used for generating a splice list, indices of where the regions have been merged, used in the region plot figure to add black seperating lines. One for each line in the bed file
	splice_full = [] # contain the length of each list in the list to compress, hence the last position of every coverage depth value list. One per row in the bed file containing regions
	for line in listToCompress:
		splice_full.append(len(line))

	splice_temp = [] # Only contains a temporary splicing value
	splice = [] # Contain lists of all splice values one list per merged region

	previous_region = [] # Previous region name 
	current_region = [] # Current region name
	combine_temp = [] # Temporary list with lines to compress
	newList=[] # Contain the new merged list from listToCompress
	
	index=-1 # Start at -1 that would be the previous element for the first element in the list, the -1 index is skipped and we move forward in the list
	regions.append('end.extra') # add one extra element to the list as for each line the previous information is stored
	listToCompress.append(['end.extra']) # add one extra elemtn to the list as for each line the previous information is stored
	

	# Merge the list to compress
	for line in regions:
		current_region = line.split('.')
		current_region = current_region[0]

		if not(previous_region ==[]):
			# As long as the previous and current region has the same name merge rows and add to the temporary combine_temp, also add the splice index to splice_temp
			if previous_region == current_region:
				combine_temp += listToCompress[index]
				splice_temp.append(splice_full[index])

			# If the Current region name is different from the previous rows region name, add the previous region to the temp list and save it in the new merged list (newlist)
			# while alsoadding the splice indices. 
			else:
				combine_temp += listToCompress[index]
				newList.append(combine_temp)
				splice_temp.append(splice_full[index])
				splice.append(splice_temp)
				# The temporary variables are emptied 
				combine_temp = []		
				splice_temp=[]		
		index+=1 # Move forward one line in the listToCompress and the splice_full list
		previous_region=current_region # Change the previous element for next iteration

	# Remove the extra element added to the lists
	regions.pop()
	listToCompress.pop()

	# Create a output list of first part of region names with no duplicates
	regions_temp=[]
	region_name=[]

	# Extract the first written region name before the '.' and save each name in the list regions temp
	for line in regions:
		l = line.split('.')
		regions_temp.append(l[0])

	#Extract duplicates from the list and keep the original order of the elements in the list
	region_name = list(OrderedDict.fromkeys(regions_temp))

	# Return the generated lists
	return region_name, newList, splice

	######################### Genrates new list of merged information rows for the validation table. The original region information list contains
	# the values chr start, stop and length. If the information comes from the same region name the values will be added horisontally 
	# list to compress=[[chr 1', 'start', 'stop'], [chr 1', 'start', 'stop'], [chr 1', 'start', 'stop'], [chr 1', 'start', 'stop']]
	# regions = ['GeneA.g.1.e', 'GeneA.g.2.e', 'GeneB.g.1.e', 'GeneC.g.1.e']

def CombineRegionInfo(listToCompress, regions):
	from collections import OrderedDict

	previous_region = []
	current_region = []
	combine_temp = []
	newList=[]
	index=-1 # Since we compare with the previous start at -1 instead of 0 to jump one step and start with the first element as a previous value

	# Since we will save the previous value on the current line iteration we need one extra made up value in both lists to make sure the last elment is saved
	regions.append('end.extra')
	listToCompress.append(['end.extra'])

	
	for line in regions:
		current_region = line.split('.')
		current_region = current_region[0]

		if not(previous_region ==[]):
			if previous_region == current_region:
				combine_temp += listToCompress[index]

			else:
				combine_temp += listToCompress[index]
				newList.append(combine_temp)
				combine_temp = []		
		index+=1
		previous_region=current_region

	# Removes the added elements to the input lists
	regions.pop()
	listToCompress.pop()

	return newList


	