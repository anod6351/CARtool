# Generate statistics table of coverage breadth values at the three given coverage depth thresholds

def stat_table(coverage_list, RegionNames, validation, phred_score, coverage_phred, X_Cut_off_list, RegionInfo, dataType):

	Tresh_1 = 0
	Tresh_2 = 0
	Tresh_3 = 0
	index=0
	s_table =[]
	s_table_phred=[]
	val_list =[]

		# Create the stat table for the coverage list without phred score filtering. 
		# The table contains the coverage breadth at three given coverage depth threshold values, X_cut_off_list.
		# If validation is turned on the validation column is added

	for line in coverage_list:
		for element in line:
			if int(element) >= int(X_Cut_off_list[0]):
				Tresh_1+=1
			if int(element) >= int(X_Cut_off_list[1]):
				Tresh_2+=1
			if int(element) >= int(X_Cut_off_list[2]):
				Tresh_3+=1


				# compute and enter the coverage depth values in the table together with region name and validation column if turned on 
		if validation:
			if float(Tresh_1)/float(len(line)) < 0.95:
				s_table.append([dataType,RegionNames[index], round(float(Tresh_1)/float(len(line)),2), round(float(Tresh_2)/float(len(line)),2), round(float(Tresh_3)/float(len(line)),2), '***'])
				val_list.append([RegionNames[index], int(round(float(Tresh_1)/float(len(line)),2)*100)]+ RegionInfo[index])
			else: 
				s_table.append([dataType,RegionNames[index], round(float(Tresh_1)/float(len(line)),2), round(float(Tresh_2)/float(len(line)),2), round(float(Tresh_3)/float(len(line)),2)])
		else:
			s_table.append([dataType, RegionNames[index], round(float(Tresh_1)/float(len(line)), 2), round(float(Tresh_2)/float(len(line)),2), round(float(Tresh_3)/float(len(line)),2)])
		Tresh_1 = 0
		Tresh_2 = 0
		Tresh_3 = 0
		index+=1


			# If phred score filtering has been turned on an additional statistics table with the filtered values is generated
			# with the same method as above.

	if phred_score:

		Tresh_1_p = 0
		Tresh_2_p = 0
		Tresh_3_p = 0
		count=0
		
		for line in coverage_phred:
			for element in line:
				if int(element) >= int(X_Cut_off_list[0]):
					Tresh_1_p+=1
				if int(element) >= int(X_Cut_off_list[1]):
					Tresh_2_p+=1
				if int(element) >= int(X_Cut_off_list[2]):
					Tresh_3_p+=1

			s_table_phred.append(['filtered',RegionNames[count], round(float(Tresh_1_p)/float(len(line)),2), round(float(Tresh_2_p)/float(len(line)),2), round(float(Tresh_3_p)/float(len(line)),2)])

			Tresh_1_p = 0
			Tresh_2_p = 0
			Tresh_3_p = 0
			count+=1

	return s_table, s_table_phred, val_list
