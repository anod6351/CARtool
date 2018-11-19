def Meanfunction(mylist, Xcov, start_pos, region_info):


	previous_element = -1 # made up start value for the previous element in mylist
	Coveragesum=0
	MeanList = []
	end_pos = 0

	count = 0 # used to caluculate end position: start + count 
	
	# If only one coverage element print the one base pair region with chromosom, start, stop, length and mean=coverage value in this case
	if len(mylist) == 1:
		MeanList.append(region_info + [int(start_pos), int(start_pos)+1, int(mylist[0])/1, 1])
		print(region_info + [int(start_pos), int(start_pos)+1, int(mylist[0])/1, 1])

	# If the list contains more than one per base coverage values, go through all elements in the list. Save the previous element each time. 
	# Divide up the regions by comparing the current and the previous elements coverage value. If both over or below keep in same group, if on different sides of the 
	# threshold coverage value save the previous and add it to the list Mean list whila starting a new region with the current element. Continue until all elements has been placed in a region in the mean list.   
	else:	
		mylist.append(0)
		for i in range(len(mylist)):

			# Make sure that the first previous element is from the actual list and not the start value of -1 
			if int(previous_element) >= 0:

				# Check if the current element coverage value is over threshold
				if int(mylist[i]) < int(Xcov):

					# If previous element also over threshold add previous element to a temporary coverage sum variable to calculate mean coverage for each new sub region
					if int(previous_element) < int(Xcov):
					 	Coveragesum += int(previous_element)
					 	count += 1
					 	# If at the end of the list save the last sub region
					 	if i == len(mylist)-1:
					 		end_pos = int(start_pos)+int(count)
					 		length = end_pos-int(start_pos)
							MeanList.append(region_info + [int(start_pos), end_pos, round(float(Coveragesum)/float(length),0), length])

					# If the previous element are instead below, save the previous region in the mean list
					else:
						Coveragesum += int(previous_element) 
						count+=1
						end_pos = int(start_pos)+int(count)
						length = end_pos-int(start_pos) 
						MeanList.append(region_info + [int(start_pos), end_pos, round(float(Coveragesum)/float(length),0), length])
						start_pos = int(end_pos)
						Coveragesum = 0
						count=0

				# Check if the current element coverage value is below threshold
				if int(mylist[i]) >= int(Xcov):

					# If the previous element is also below threshold add previous value to the coverage sum
					if int(previous_element) >= int(Xcov):
						Coveragesum += int(previous_element)
						count+=1
						# If the end of the list is reached save the last region
						if i == len(mylist)-1:
							end_pos = int(start_pos)+int(count)
							length = end_pos-int(start_pos)
							MeanList.append(region_info + [int(start_pos), end_pos, round(float(Coveragesum)/float(length),0), length])

					# If the previous and current is on opposite sides of the threshold save the new sub region and start a new region with the current element
					else:
						Coveragesum += int(previous_element)
						count+=1
						end_pos = int(start_pos)+int(count)
						length = end_pos-int(start_pos)
						MeanList.append(region_info + [int(start_pos), end_pos, round(float(Coveragesum)/float(length),0), length])
						start_pos = int(end_pos)
						Coveragesum = 0
						count=0

			previous_element = int(mylist[i]) # Move the previous element forward for the next iteration
			i+=1 # index for the coverage value list (mylist)

		mylist.pop()

	return MeanList






