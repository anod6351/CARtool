def StrandDifference(positive_table, negative_table):

	# Declare new table variables, that will be filled with region name, three coverage breadth values and a difference column between the two strand specific tables. 
	positive_final=[]
	negative_final=[]
	
	# contain the difference in percentage between the first coverage breadth values in each strand specific table
	diff = 0

	for pos_line in positive_table:
		for neg_line in negative_table:
			diff = abs(int(neg_line[2])-int(pos_line[2]))
			positive_final.append(pos_line[1:5] + [diff])
			negative_final.append(neg_line[1:5] + [diff])

	return positive_final, negative_final

