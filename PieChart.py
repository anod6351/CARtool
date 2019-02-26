
# Generate pie chart

def Generate_Pie(Detailed_list, Threshold, index, fig, title, textSize):

	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	import tkinter

	# Count values below and over threshold 
	count_under = 0
	count_over = 0

	for element in Detailed_list:
		if int(element) < int(Threshold):
			count_under+=1
		else:
			count_over+=1

	# Generate a pie chart with the count over and under threshold. Blue for over and yellow for the values bellow.
	ax = fig.add_subplot(2,2,index)
	ax.pie([count_over,count_under], colors=['blue', 'orange']) 
	ax.axis('equal')
	ax.set_title(str(title), fontsize=textSize)
	ax.legend([str(count_over) +' bp. > ' + str(int(Threshold)) +'X', str(count_under) + ' bp. < '+ str(int(Threshold)) + 'X'], loc='lower center', bbox_to_anchor=(0.5, 0), fontsize='xx-small')

	return fig
