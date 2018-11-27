# Program lancher for Coverage Analysis Report, CAR

import argparse
import sys
import pybedtools
import csv
import os

##############################  Parameter input #########################################

parser = argparse.ArgumentParser()

# Mandatory input, creates:
# 1. A coverage depth short list with subregions below or above the first coverage depth threshold in the X_cut_off_list  
# 2. A coverage breadth table with percentage at or over the values in the X_cut_off_list
# 3. A Logg file
parser.add_argument("-a","--Regions", help="Write the path to the Bed file")
parser.add_argument("-b","--Reads", help="Write the path to the Bam file")
parser.add_argument("-c","--X_Cut_off_list", nargs='+', help="Write coverage depths thresholds, the first value will be used for coverage depth analytics and all three for coverage breadth")

# Additional optional input settings
parser.add_argument("-p", "--phred_score", nargs='+', help="Choose a phred score threshold for additional row in the statististics table")
parser.add_argument("-v", "--validation", action="store_true", help="turn validation list option on")
parser.add_argument("-k", "--combineRows", action="store_true", help="turn per gene option on, combines rows in the bedfile from the same gene for example")
parser.add_argument("-s","--strandSpecific", action ="store_true", help="Create an additional statistics table with forward and reverse read coverage calculated seperatly")
parser.add_argument("-f", "--figures", action="store_true", help="Create figures") 
#parser.add_argument("-t", "--hotspot", help="Create figures")
parser.add_argument("-t", "--hotspot", action='store_true', help="Create figures")
parser.add_argument("-e", "--executedBy", nargs='+', help="Enter name of person that run the coverage tool")
#parser.add_argument("-l", "--lowRegions", help="List of known to be low regions")
parser.add_argument("-l", "--lowRegions", action='store_true', help="List of known to be low regions")
parser.add_argument("-m", "--ExonTranscript", action="store_true", help="adds exon and transcript information to the low coverage mean region list")
parser.add_argument("-d", "--detailedCoverage", action="store_true", help="Outputs a per position coverage depth list of the region bed file")
args = parser.parse_args()

ProgramVersion = "Program Verion 1.0"

# Low region and hotspot list for testing during construction: 
lowRegions = []
hotspots = []
if args.lowRegions:
	lowRegions = [['CSF3R', 'chr1', 36931697, 36931787], ['CSF3R', 'chr1', 36932117, 36932432], ['NRAS', 'chr1', 115256540, 115256599]]
if args.hotspot:
	hotspots=[['CSF3R', 'chr1', 36931697, 36931698], ['CSF3R', 'chr1', 36931786, 36931787],  ['NRAS', 'chr1', 115256580, 115256581]]


########################### Check that all mandatory parameter input is given

if args.Regions == None or args.Reads == None or len(args.X_Cut_off_list) < 3:
	sys.exit("Not enough input: example --Regions myRegions.bed --Reads myReads.bam --Xcov 20 --X_Cut_off_list 10 20 30. Use help for more information")

if not(args.executedBy):
	sys.exit("Please enter your name, the name of who run the program will be stored in the logg file")

############################# RUN THE PROGRAM ###################################################

print("Running coverage analysis ...")


# Create the an output folder for the analytic report
#command = "mkdir CAR_output"
#os.popen(command)

####################### Create per base pair coverage list ###################################################################
# The script uses bedtools and formats the output so that each row contains all coverage values from the corresponding region (row) in the bedfile.  

if args.phred_score==None or not(str(args.phred_score[0])=='all'):


	# Deal with the region file being one based at the start position
	regions_file = open(args.Regions, 'r')
	regions_fixed_file = [] 
	for line in regions_file:
		element=line.split('\t')
		regions_fixed_file.append([element[0], int(element[1])-1, element[2], element[3]])
	Region_fixed = pybedtools.BedTool(regions_fixed_file)


	# Create bedtool objects
	print("Creating BedTool objects ...") 

	Reads = pybedtools.BedTool(args.Reads)
	Regions =pybedtools.BedTool(args.Regions)

	# Set the data type variable to raw data, this will be shown in the statistics table
	dataType = 'Raw'

	# Create per base/amplicon coverage
	print("Creating per base/amplicon coverage ...")
	import Detailed_generator
	detailed_list = Detailed_generator.detail(Region_fixed, Reads) ############ Regions_fixed if 1 based start position


########################## Phred score filtering of bam file (OPTIONAL) ################################################################################### 

# Use samtools to filter the bam file with the given threshold value. This can either be used for the whole analysis or only as an additional field in the statistics table. 
# With the input -p all 20, the bam file with reads will be filtered and used for all calclations. While -p 20 will filter the bam file and add to the statistics table with the non filtered bam file as mainly used. 

detailed_list_phredfilter=[]
phred_score = False

# Check the type of phred score setting, save the phred score as a variable
if args.phred_score:
	print("Creates a phred score filtered bam file for the statistics table ...")
	if not(str(args.phred_score[0]) == 'all'):
		phred_score = args.phred_score[0]
	else:
		phred_score=args.phred_score[1]
		dataType= 'Phred filtered' # this will be shown in the statistics table

		#### Deal with the start position index being 1 based in the bed file...

		regions_file = open(args.Regions, 'r')
		regions_fixed_file = [] 

		for line in regions_file:
			element=line.split('\t')
			regions_fixed_file.append([element[0], int(element[1])-1, element[2], element[3]])

		Region_fixed = pybedtools.BedTool(regions_fixed_file)

	command = "samtools view -bq " + str(phred_score) + " " + str(args.Reads) + " > CAR_output/filtered_bam.bam"
	p = os.popen(command)
	p.wait()
	
	Reads_filtered = pybedtools.BedTool(open("CAR_output/filtered_bam.bam", "r"))

	# Compute per base coverage lists

	if str(args.phred_score[0])=="all":
		Regions =pybedtools.BedTool(args.Regions)
		import Detailed_generator
		detailed_list = Detailed_generator.detail(Region_fixed, Reads_filtered) #######Depends on the bed file
	else:
		detailed_list_phredfilter = Detailed_generator.detail(Region_fixed, Reads_filtered) #######Depends on the bed file

############################# Filter bam file with reads to only contain + or - strand (OPTIONAL) #################################################

# These strand specific detailed coverage lists will be used to create statistic tables containing coverage breadth values at the three threshold values
if args.strandSpecific:
	print("Creates strand specific bam files for the additional statistics table ...")
	command_1 = "samtools view -F 0x10 -b " + str(args.Reads) + " > CAR_output/positive_strand.bam"
	command_2 = "samtools view -f 0x10 -b " + str(args.Reads) + " > CAR_output/negative_strand.bam"
	os.popen(command_1)
	os.popen(command_2)

	# Create per position coverage for the strand specific bam files 
	Reads_positive = pybedtools.BedTool(open("CAR_output/positive_strand.bam", "r"))
	detailed_list_positive = Detailed_generator.detail(Regions, Reads_positive)
	Reads_negative = pybedtools.BedTool(open("CAR_output/negative_strand.bam", "r"))
	detailed_list_negative = Detailed_generator.detail(Regions, Reads_negative)




###################### CombineRows (OPTIONAL) ################################################################################################## 
# Merge rows in the per position coverage lists (detailed lists) so that ALL coverage values from 
# the same, for example, gene will be in the same row. This is essential for per gene calculations. The output will be a new formated 
# detailed coverage result list and a combined region name list. The region name in the bedfile must be seperated with a '.' .Example: gene1.exon.2 etc. 
# for this function

det_list_phred_formated =[]
splice=[] # Gives the indices of the subregions within the combined rows
Info=[] # Gives the chromosome, start and stop for each subregion in the combined row region

# Create a list that only containing the region name column
Regions_list=[]
for line in Regions:
	Regions_list.append(line[3])

# Combine rows in the bedfile that has the same region name before the first dot as separator. Example: the rows Gene1.Exon.2 and Gene1.Exon.3 will be merged
if args.combineRows:
	import CombineRows_generator
	print("Combine rows ...")

	# Formats the detailed coverage list by combining rows
	RegionNames, detailed_list_formated, splice = CombineRows_generator.CombineRowsList(detailed_list, Regions_list)

	# If phred score filtering option is turned on the phred filtered detailed coverage list is fomrated by merging rows  
	if args.phred_score:
		PhredRegionNames, det_list_phred_formated, splice_phred = CombineRows_generator.CombineRowsList(detailed_list_phredfilter, Regions_list)

	if args.strandSpecific:
		# If strand specific option turned on  the negative and positive strand coverage detail lists are formated by merging rows
		PosRegionNames, detailed_list_positive_formated, splice_pos = CombineRows_generator.CombineRowsList(detailed_list_positive, Regions_list)
		NegRegionNames, detailed_list_negative_formated, splice_neg = CombineRows_generator.CombineRowsList(detailed_list_negative, Regions_list)



###################### Create mean coverage regions list ##############################################################################################
# The sub regions in the list are all under or over the first coverage depth threshold. In the Mean coverage short list only the sub regions below the 
# chosen coverege depth threshold are kept. 
# Example: detailed coverage list = [1, 1, 1, 4, 5, 2, 2] and if the coverage threshold is = 3 
# 	=> Mean coverage list = [[start = 0, stop = 3, mean = 1], [3, 5, 4.5], [5, 7, 2]]
# 	=> The short mean list = [[0,3,1],[5,7,2]] 

print("Generating Mean coverage short lists ...")


# create a list with only start positions from each line in the bedfile,
# this list is used to compute new start and stop positions for the new subregions
Start_pos_list =[]
for line in Regions:
	Start_pos_list.append(int(line[1]))

## Add columns to the mean region short list
info_list =[] 
info_temp=[]
mean_index=0

if args.ExonTranscript: # Adds exon number, transcript and chromosome
	for element in Regions_list:
		info_temp = element.split('.')
		info_temp = [info_temp[2], info_temp[4], info_temp[5]]
		info_list.append(info_temp)
	mean_index=5

else: # Only add the chromosome not exon and transcript
	for element in Regions_list:
		info_temp = element.split('.')
		info_temp = [info_temp[1]] ##################5 innan, annars 1 
		info_list.append(info_temp)	
	mean_index=3


# Create the Mean coverage list
import Meanlist_generator
Row_temp = []
MeanCov_List = []
index = 0

# Compute the mean sub regions for each row in the detailed list and create the Mean Coverage list, this is used to create the mean coverage region short list  
for line in detailed_list:
	Row_temp = Meanlist_generator.Meanfunction(line, args.X_Cut_off_list[0], Start_pos_list[index], info_list[index])
	MeanCov_List.append(Row_temp)
	index += 1

# Merge the rows in the mean coverage list generated above, only if the combine rows option is activated
if args.combineRows:
	import CombineRows_generator
	MeanNames, MeanCoverage_formated, splice_meanList = CombineRows_generator.CombineRowsList(MeanCov_List, Regions_list) 


##### Create a short list of mean coverage regions under the coverage depth threshold

MeanCov_shortList=[]
MeanCov_shortList_temp = []
MeanCov_shortList_sublist_temp=[]
Mean_Short_input_list = []
Mean_Short_w_sublists = []

# Decide on what Mean coverage region list to use for the short list, depends on if combine rows is activated or not 
if args.combineRows:
	Mean_full_list = MeanCoverage_formated
else:
	Mean_full_list = MeanCov_List
	MeanNames = Regions_list

# Generate the short mean subregion list by extracting low regions from the Mean coverage region list
for line in Mean_full_list:
	for element in line:
		if int(element[mean_index]) < int(args.X_Cut_off_list[0]): # Check if the subregion has a mean below the coverage depth threshold. Only the "bad" coverage subregions are stored in the short list
			MeanCov_shortList_temp+=(element)
			MeanCov_shortList_sublist_temp.append(element)

	MeanCov_shortList.append(MeanCov_shortList_temp)
	Mean_Short_w_sublists.append(MeanCov_shortList_sublist_temp)
	MeanCov_shortList_temp = []
	MeanCov_shortList_sublist_temp=[]


################################# Known low regions as input (OPTIONAL)
if args.lowRegions:
	name_lowRegions=[]

	for line in lowRegions:
		name_lowRegions.append(line[0])

	index=0
	count=0
	index_low=0
	LowRegionIndex=[]


	for name in MeanNames:

		for low_name in name_lowRegions:
			if str(name) == str(low_name):

				Low_region_compare = lowRegions[index_low]
				index_low+=1

				for line in Mean_Short_w_sublists[index]:

					if str(line[0]) == Low_region_compare[1] and line[mean_index-2]==Low_region_compare[2] and line[mean_index-1]==Low_region_compare[3]:

						LowRegionIndex.append([index, count])
					count+=1

				count=0	
		index+=1

	Mean_Short_with_lowRegionInfo = []
	Temp_row=[]
	index_line= 0
	index_low=0
	index_element = 0

	for line in Mean_Short_w_sublists:
		for element in line:
			if index_low < len(LowRegionIndex):
				Low_line = LowRegionIndex[index_low]
				if int(index_line) == int(Low_line[0]) and int(index_element) == int(Low_line[1]):
					Temp_row += element + ['Yes']
					index_low+=1
				else:
					Temp_row += element + ['No']
			else:
				Temp_row += element + ['No']

			index_element+=1

		Mean_Short_with_lowRegionInfo.append(Temp_row)
		Temp_row=[]
		index_line+=1
		index_element=0

################################ Hotspot list as input (OPTIONAL)

if args.hotspot:
	name_hotspots=[]

	for line in hotspots:
		name_hotspots.append(line[0])

	index=0
	count=0
	index_hotspot=0
	hotspotIndex=[]
	counter_start_pos=0


	for mean_name in MeanNames:

		for hot_name in name_hotspots:
			if str(mean_name) == str(hot_name):

				hotspots_compare = hotspots[index_hotspot]
				index_hotspot+=1
				for line in Mean_Short_w_sublists[index]:
					if str(line[0]) == hotspots_compare[1] and line[mean_index-2] <= hotspots_compare[2] and line[mean_index-1] >= hotspots_compare[3]:
						# Save the row index of the name of the region with a hotspot together with the sub region index and the start position of the hotspot
						counter_start_pos = int(hotspots_compare[2]) - int(line[mean_index-2])
						print(counter_start_pos)
						hotspotIndex.append([index, count, counter_start_pos])
					count+=1

				count=0	
		index+=1

			 

###################### Create statistics table #####################################################################################################
# Uses the detailed coverage list (formated if combineRows setting is turned on) reports coverage breadth over the choosen coverage depth cut off values. 
#The statistics table can also contain results from a phred filtered bam file and an optional validation column. 
# The validation option checks if 95% coverage bredth or more at the first X cut off value. If under 95% the column is marked with **** 
# and the region is added to the validation list. If the strand specific option is turned on an additional statistics table with strand specific reads is generated.

print("Generating Statistics table ...")

# If validation option is turned on, create a list with region info for the validation table. With the columns chr, start, stop, length
Region_info=[]

if args.validation:
	for line in Regions:
		Region_info.append([line[0], line[1], line[2], int(line[2])-int(line[1])])

	# If combine rows activated, merge rows in the regon info list
	if args.combineRows:
		import CombRegionInfo
		Region_info = CombRegionInfo.CombineRegionInfo(Region_info, Regions_list)

# Decide on what coverage lists to use for the statistics table
if args.combineRows:
	detailed_list_stat = detailed_list_formated
	detailed_phred_stat = det_list_phred_formated

	if args.strandSpecific:
		detailed_positive = detailed_list_positive_formated
		detailed_negative = detailed_list_negative_formated
	Stat_table_names = RegionNames
else:
	detailed_list_stat = detailed_list
	detailed_phred_stat = detailed_list_phredfilter

	if args.strandSpecific:
		detailed_positive = detailed_list_positive
		detailed_negative = detailed_list_negative
	Stat_table_names = Regions_list

	
import Statistics_generator
# Compute statistics table 	
stat_table, stat_table_phred, validation_list = Statistics_generator.stat_table(detailed_list_stat, Stat_table_names, args.validation, phred_score, detailed_phred_stat, args.X_Cut_off_list, Region_info, dataType)    
# Compute the strand specific statistics tables
if args.strandSpecific:
	stat_table_positive, stat_table_phred_positive, validation_list_positive = Statistics_generator.stat_table(detailed_positive, Stat_table_names, args.validation, False, [], args.X_Cut_off_list, Region_info, dataType)
	stat_table_negative, stat_table_phred_negative, validation_list_negative = Statistics_generator.stat_table(detailed_negative, Stat_table_names, args.validation, False, [], args.X_Cut_off_list, Region_info, dataType)

	# Add a difference column between the first column coverage breadth values in the positive and negative statistics table
	import StrandSpecific_Diff
	Stat_table_positive_final, Stat_table_negative_final = StrandSpecific_Diff.StrandDifference(stat_table_positive, stat_table_negative)


################### Create the logg ###############################################################################################################

print("Creates the logg file ...")
Logg=[]
import datetime

# Program version, variable created at the top of the script:
Logg.append([ProgramVersion])

# Add date the program was run
Logg.append(["Date: ", str(datetime.datetime.now())])

# Add the name of the person running the program 
Logg.append(["Coverage analysis run by: ", args.executedBy])

# Mandatory input:
# Add the region and reads files used
Region_fileName = args.Regions
Reads_fileName = args.Reads
Logg.append(["Region file (bedfile): ", Region_fileName])
Logg.append(["Read file (bamfile): ", Reads_fileName])
# Add the chosen coverage depth thresholds
Logg.append(["Coverage depth thresholds: ", args.X_Cut_off_list[0], args.X_Cut_off_list[1], args.X_Cut_off_list[2]])

# Optional input:

# Phres score
if args.phred_score:
	if str(args.phred_score[0])=="all":
		Logg.append(["Choosen phred score: ", args.phred_score[1], "Used in all calculations"])
	else:
		Logg.append(["Choosen phred score: ", args.phred_score[0], "Used only in an additional row in the statistics table"])

# validation option
if args.validation:
	Logg.append(["Validation list generated"])

# Combine rows
if args.combineRows:
	Logg.append(["Combined region rows in bedfile"])

# Generation of figures
if args.figures:
	Logg.append(["Figures generated"])
	# If hotspots been added to the region figure
	if args.hotspot: 
		Logg.append(["Hotspots added to the region figure"]) # add what hotspot list that was used

# If exon and transcript has been added to the mean coverage list
if args.ExonTranscript:
	Logg.append(["Exon and transcript added to mean region lists"])

# If a list of low regions was given and an additional column was added to the mean coverage list
if args.lowRegions:
	Logg.append(["Low region list added", args.lowRegions])

# If the detailed list of coverage depth per position was saved as output
if args.detailedCoverage:
	Logg.append(["Per base coverage depth list generated"])

Logg.append(["Calculations:"])

#Calulate the total mean coverage and add to the logg file
Total_Mean_Cov = 0
count=0

for line in Mean_full_list: 
	for element in line:
		Total_Mean_Cov += element[mean_index]
		count+=1

Total_Mean_Cov = round(float(Total_Mean_Cov)/float(count),2)
Logg.append(["Mean coverage width:", str(Total_Mean_Cov) + ' X'])
	
# Calculate the mean coverage breadth for the columns in the statistics table and add to the logg file
stat_temp1 = 0.0
stat_temp2 = 0.0
stat_temp3 = 0.0

for line in stat_table:
	stat_temp1 += line[2]
	stat_temp2 += line[3]
	stat_temp3 += line[4]
Logg.append(["Mean coverage breadth at: ",  args.X_Cut_off_list[0], args.X_Cut_off_list[1], args.X_Cut_off_list[2]])
Logg.append(["", round(float(stat_temp1)/float(len(stat_table)),2),round(float(stat_temp2)/float(len(stat_table)),2), round(float(stat_temp3)/float(len(stat_table)),2)])


######################################## Figures ######################################################################################################
if args.figures == True:

	# Diside what detailed coverage list and region name list to use, depending on if combined rows is turned on or not 
	if args.combineRows:
		list_for_figures = detailed_list_formated
		Names_for_figures = RegionNames
	else:
		list_for_figures = detailed_list
		Names_for_figures = Regions_list

	##### PIE CHART #######################################################################################################################################
	print("Generating figures ...")

	# One pie chart per row in bedfile (region file) or one per combined region. The pie charts are only created for regions that has some subpart below the coverage depth threshold, X_cut_off_list[0]
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	import Tkinter
	import PieChart

	fig = plt.figure()
	index_order = 1 # Used to place the figure in the 2*2 grid in the figure and as an indicator pie charts per PDF 
	image_count = 1 # Used in the figure names as an index
	cov_under = 0 # Used to check if any elements in the each row in the detailed coverage list has values below the threshold, X_cut_off[0] 
	Name_index = 0 # Used for selecting the right name for the pie charts that are generated

	if args.combineRows:
		text_size = 'small'
	else:
		text_size = 'xx-small'

	# Create the pie charts
	for line in list_for_figures:

		# Check if any elements in the each row in the detailed coverage list has values below the threshold, X_cut_off[0] 
		for element in line:
			if int(element) < int(args.X_Cut_off_list[0]):
				cov_under+=1

		# Only create pie charts with lines that has a some position under the coverage threshold
		if not(int(cov_under) == 0):	
			fig = PieChart.Generate_Pie(line, args.X_Cut_off_list[0], index_order, fig, Names_for_figures[Name_index], text_size)

			# when 4 pie charts exists in the figure they are saved in an indexed pdf file
			if index_order==4:
				# Save the pie charts as a PDF 
				fig.savefig('CAR_output/Cov_pie_' + str(image_count) + '.pdf')
				fig = plt.figure()
				index_order=0
				image_count+=1
			index_order+=1
		cov_under=0
		Name_index+=1

	# Make sure that the last image is printed even if there is less that 4 pie charts in it
	if not(fig == None): 
		fig.savefig('CAR_output/Cov_pie_' + str(image_count) + '.pdf')

############################################# Bar plot

	if args.combineRows:

		# Create a list with subpart names that is used as x-axis labels for each bar and a .....
		import Subpart_names 
		subpartNames, detailed_list_formated_bar = Subpart_names.bar_names_generator(detailed_list, Regions_list)
		
		# Create the bar plot, bar plots are only generated if the combined regions contain some element with a coverage value below the threshold, x_cut_off_list[0]
		import bar_plot
		fig = plt.figure()
		index_order = 1
		Name_index = 0
		image_count = 1
		bellow = 0

		for line in detailed_list_formated_bar:
			for sublist in line:
				for element in sublist:
					if int(element) < int(args.X_Cut_off_list[0]):
						bellow+=1

			if int(bellow) > 0:
				fig = bar_plot.bar_plot_generator(line, args.X_Cut_off_list[0], index_order, fig, subpartNames[Name_index], Names_for_figures[Name_index])
				
				# When 4 bar diagrams exists in fig they are saved in a indexed pdf file
				if index_order==4:
					fig.savefig('CAR_output/Cov_bar_' + str(image_count) + '.pdf')
					fig = plt.figure() # Create a new figure to draw the pie charts in
					index_order=0
					image_count+=1
				index_order+=1
			Name_index+=1
			#print(bellow)
			bellow=0

		# Make sure that the last image is printed even if there is less that 4 pie charts in it
		if not(fig == None): 
			fig.savefig('CAR_output/Cov_bar_' + str(image_count) + '.pdf')

	###### Per region coverage plot ###############################################################################################################

	print("Generating region figure ...")


	############## Create text list for regions, this will be printed next to the region plots
	Text_list_temp = []
	Text_list = []
	Figure_text_full =[]

	if args.combineRows:
		Figure_text_full = MeanCoverage_formated
	else:
		Figure_text_full = MeanCov_List

	for line in Figure_text_full:
		for element in line:
			if float(element[mean_index]) < float(args.X_Cut_off_list[0]):
				Text_list_temp.append(element)
		Text_list.append(Text_list_temp)
		Text_list_temp = []


	########################################################################################

	import Region_position_plot
	Name_index = 0
	index_order = 1
	image_count = 1
	count_under = 0
	#count=1
	only_under_index = 0
	data = []
	column =[]
	rows=[]
	fig = plt.figure()
	splice_index=0
	hotspots_arrows=[]
	hotspot_index=0

	if splice == None: 
		splice=[]

	for line in list_for_figures:
		for element in line: 
			if int(element) < int(args.X_Cut_off_list[0]):
				count_under+=1

		if int(count_under) > 0:

			if args.hotspot:
				if hotspot_index < len(hotspots):
					for hs_line in hotspotIndex:
						if int(Name_index) == int(hs_line[0]):
							hotspots_arrows.append(hs_line[2])
					hotspot_index+=1

			if args.combineRows:
				fig = Region_position_plot.Region_generator_plot(line, Names_for_figures[Name_index], args.X_Cut_off_list[0], index_order, fig, splice[splice_index], hotspots_arrows)
			else:
				fig = Region_position_plot.Region_generator_plot(line, Names_for_figures[Name_index], args.X_Cut_off_list[0], index_order, fig, False, hotspots_arrows)
			hotspots_arrows=[]
			index_order+=1

			ax = fig.add_subplot(2,2,index_order)

			for element in Text_list[only_under_index]:
				data.append(element)

			if not(args.ExonTranscript):	
				column = ['Chr','start', 'stop', 'mean','length']
			else:
				column = ['Chr', 'exon', 'transcript','start', 'stop', 'mean', 'length']

			for i in range(len(data)):
				rows.append('R' +  str(i+1))

			# If low regions are sent as input, add red color to any line in the table containg a known to be low region
			color=[]
			line_index=0
			index_low=0
			low_line=[]
			if args.lowRegions:
				for row in data:
					if int(index_low) < int(len(LowRegionIndex)):
						Low_line=LowRegionIndex[index_low]
						if int(line_index) == int(Low_line[1]):
							color.append('tomato')
							index_low+=1
						else:
							color.append('white')
					else:
						color.append('white')
					line_index+=1

			if color==[]:
				fig_table = ax.table(cellText= data, rowLabels = rows, colLabels = column, loc='center')
			else:
				fig_table = ax.table(cellText= data, rowLabels = rows, colLabels = column, loc='center', rowColours=color)
			data=[]
			rows=[]
			ax.axis('off')

					
			if index_order == 4:
				fig.savefig('CAR_output/RegionCoverage_' + str(image_count) + '.pdf')
				fig = plt.figure()
				index_order = 0
				image_count+=1

			index_order+=1

		Name_index+=1
		splice_index+=1
		count_under = 0
		only_under_index+=1


	#Make sure that the last region plot image is printed even if there is less than 4 region plots in it
	if not(fig == None): 
		fig.savefig('CAR_output/RegionCoverage_' + str(image_count) + '.pdf')

	

################## Save list and tables to csv files ###################################################################################################

print("Saving to csv files ...")

###### Save the non strand specific statistics table to a csv file

index=0
with open('CAR_output/Stat_table.csv', 'w') as myfile:
	wr = csv.writer(myfile)

	if args.validation:
		wr.writerow(['Data type','Region name', str(args.X_Cut_off_list[0]) +'X', str(args.X_Cut_off_list[1])+'X', str(args.X_Cut_off_list[2])+'X', 'Validation'])
	else:
		wr.writerow(['Data type','Region name', str(args.X_Cut_off_list[0]) +'X', str(args.X_Cut_off_list[1])+'X', str(args.X_Cut_off_list[2])+'X'])
	for row in stat_table:
		wr.writerow(row)
		if args.phred_score and not(str(args.phred_score[0])) == "all":
			print(index)
			wr.writerow(stat_table_phred[index])
			index+=1

###### Save the positive strand specific statistics table to a csv file
if args.strandSpecific:
	index=0
	with open('CAR_output/Stat_table_positiveStrand.csv', 'w') as myfile:
		wr = csv.writer(myfile)
		wr.writerow(['Region name', str(args.X_Cut_off_list[0]) +'X', str(args.X_Cut_off_list[1])+'X', str(args.X_Cut_off_list[2])+'X', 'strand difference'])
		for row in Stat_table_positive_final:
			wr.writerow(row)

###### Save the negative strand specific statistics table to a csv file

	index=0
	with open('CAR_output/Stat_table_negativeStrand.csv', 'w') as myfile:
		wr = csv.writer(myfile)
		wr.writerow(['Region name', str(args.X_Cut_off_list[0]) +'X', str(args.X_Cut_off_list[1])+'X', str(args.X_Cut_off_list[2])+'X', 'strand difference'])
		for row in Stat_table_negative_final:
			wr.writerow(row)


###### save Mean coverage short list to csv file

index=0
with open('CAR_output/MeanCoverageShortList.csv', 'w') as myfile:
	wr =csv.writer(myfile)
	if args.ExonTranscript:
		if args.lowRegions:
			wr.writerow(['Region name', 'chr', 'exon', 'transcript', 'start', 'stop', 'Mean Coverage', 'length', 'known'])
		else:
			wr.writerow(['Region name', 'chr', 'exon', 'transcript', 'start', 'stop', 'Mean Coverage', 'length'])
	else:
		if args.lowRegions:
			wr.writerow(['Region name', 'chr', 'start', 'stop', 'Mean Coverage', 'length', 'known'])
		else: 
			wr.writerow(['Region name', 'chr', 'start', 'stop', 'Mean Coverage', 'length'])

	for row in MeanCov_shortList:
		if not(row == []):
			if args.lowRegions:
				line = [MeanNames[index]] + Mean_Short_with_lowRegionInfo[index]
				wr.writerow(line)
			else:
				wr.writerow([MeanNames[index]] +row)

		index+=1

###### save Mean Coverage Full List to csv file

index=0
with open('CAR_output/MeanCoverageFullList.csv', 'w') as myfile:
	wr= csv.writer(myfile)

	if args.ExonTranscript:
		if args.lowRegions:
			wr.writerow(['Region name', 'chr', 'exon', 'transcript', 'start', 'stop', 'Mean Coverage', 'length', 'known'])
		else:
			wr.writerow(['Region name', 'chr', 'exon', 'transcript', 'start', 'stop', 'Mean Coverage', 'length'])
	else:
		if args.lowRegions:
			wr.writerow(['Region name', 'chr', 'start', 'stop', 'Mean Coverage', 'length', 'known'])
		else: 
			wr.writerow(['Region name', 'chr', 'start', 'stop', 'Mean Coverage', 'length'])

	for row in MeanCov_List:
		row = [Regions_list[index]] + row
		wr.writerow(row)
		index+=1

##### Save the logg file to csv file

with open('CAR_output/Logg.csv', 'w') as myfile:
	wr = csv.writer(myfile)
	wr.writerow(['LOGG'])
	for line in Logg:
		wr.writerow(line)

# Save the validation list for the original statistics table to a csv file

if args.validation:

	with open('CAR_output/Validation_list.csv', 'w') as myfile:
		wr = csv.writer(myfile)
		wr.writerow(['Region name', '% Coverage at: ' + str(args.X_Cut_off_list[0])+ ' X', 'Chr', 'Start', 'Stop', 'length']) 
		for line in validation_list:
			wr.writerow(line)

# Saves the per base coverage depth list to a csv file

if args.detailedCoverage:
	index=0
	Region_info_det=[]
	for line in Regions:
		Region_info_det.append([line[3], line[0], line[1], line[2]])

	with open('CAR_output/PerBaseCoverage.csv', 'w') as myfile:
		wr =csv.writer(myfile)
		wr.writerow(['Region name', 'Chr', 'Start', 'Stop','Coverage Depth'])
		for line in detailed_list: 
			line = Region_info_det[index] + line
			wr.writerow(line)
			index+=1
