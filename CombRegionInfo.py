
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


	