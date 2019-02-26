# Generate the postive strand and negative strand statistics table of coverage breadth values at the three given coverage depth threshold values, along with a column defining the difference between the two tables for the first coverage depth threshold.
def StrandDifference(positive_table, negative_table):

	# Declare new table variables, that will be filled with region name, three coverage breadth values and a difference column between the two strand specific tables. 
	positive_final=[]
	negative_final=[]
	
	# contain the difference between the first coverage breadth values in each strand specific table
	diff = 0
	pos_index=0
	pos_line=[]
	
	for neg_line in negative_table:

		pos_line = positive_table[pos_index]
		diff = round(abs(float(neg_line[2])-float(pos_line[2])),2)

		positive_final.append(pos_line[1:5] + [diff])
		negative_final.append(neg_line[1:5] + [diff])

		pos_index+=1

	return positive_final, negative_final
