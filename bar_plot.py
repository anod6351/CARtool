# bar plot 
def bar_plot_generator(data_list, C_threshold, index, fig, RegionNames, titleName):
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt 

	ax = fig.add_subplot(2,2,index)
	data=[]

	# Checks if the exon numbers are added in growing order else reverse the data list with each regions covareg values and the region names  

	#if isinstance(RegionNames[0], int) and isinstance(RegionNames[1], int):
	#	print('yes')
	if int(len(RegionNames)) > 1 and int(RegionNames[1]) < int(RegionNames[0]):
			data=list(reversed(data_list))
			RegionNames=list(reversed(RegionNames))
	else:
		data=data_list

	
	width = 0.35 # width of each position representation in the graph
	over_count=0 # variable to store the temporary number of base pairs above the threshold
	under_count=0 # variable to store the temporary number of base pairs below the threshold
	data_over = [] # list of amounts over threshold values
	data_under = [] # list of amounts under threshold values

	# Compute the fraction of values under and over coverage depth threshold for all regions in the data list, (one line in the bedfile or several merged lines if combine rows was activated)
	for line in data:
		for element in line:
			if int(element) < int(C_threshold):
				under_count+=1
			else:
				over_count+=1
		# Only adds the bar with a score below 
		#if int(under_count) > 0:
		data_over.append(round((float(over_count)/float(len(line))),1) * 100)
		data_under.append(round(float(under_count)/float(len(line)),1) * 100)
		# Reset temporary variables to 0 for next line calculations
		under_count=0
		over_count=0

	# Creat list with x coordinate positions for each bar
	x_position =[] 
	for i in range(len(data_over)):
		x_position.append(i+1)

	#Add the fraction data (amount of below/above values), x coordinates for the bars, size and color of each part of the bars
	ax.bar(x_position, data_over, width, color='blue')
	ax.bar(x_position, data_under, width, color='orange', bottom=data_over)

	# Scale of X axis, if the number of bars in the plot is 10 or below the x axis is between 0 and 10. If the bars are more than 10 the x axis is scaled accordingly  
	if len(data_list)>10:
		ax.set_xlim(0,len(data_over)+1)
	else:
		ax.set_xlim(0,10)

	
	ax.set_ylim(0,100) # Set the limit on the y axis to 100, 100% is always the maximum
	ax.set_title(titleName) # Set the name of the figure, placed above 
	if len(data_list)> 10:
		plt.xticks(x_position, RegionNames, fontsize='xx-small')
	else:
		plt.xticks(x_position, RegionNames) # Add the region names or exon number for example to the x axis

	return fig