
# Generate pie chart

def Generate_Pie(Detailed_list, Threshold, index, fig, title, textSize):

	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	import Tkinter

	count_under = 0
	count_over = 0
	print(Detailed_list)
	for element in Detailed_list:
		if int(element) < int(Threshold):
			count_under+=1
		else:
			count_over+=1

	print(count_over)
	print(count_under)
	ax = fig.add_subplot(2,2,index)
	ax.pie([count_over,count_under], colors=['blue', 'orange']) 
	ax.axis('equal')
	ax.set_title(str(title), fontsize=textSize)
	ax.legend([str(count_over) +' bp. > ' + str(int(Threshold)) +'X', str(count_under) + ' pos. < '+ str(int(Threshold)) + 'X'], loc='lower center', bbox_to_anchor=(0.5, 0), fontsize='xx-small')

	return fig

