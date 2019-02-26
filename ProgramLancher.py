# Program lancher for Coverage Analysis Report, CAR

import argparse
import csv
import subprocess

##############################  Parameter input #########################################

parser = argparse.ArgumentParser()

# Mandatory input, creates:
# 1. Two coverage depth lists (short and full lists) of subregions below or above the first coverage depth threshold in the X_cut_off_list  
# 2. A coverage breadth table with percentage of values equal and above the thresholds in the X_cut_off_list
# 3. A Logg file, the name of who run CAR is needed as input
# 4. Create a folder where the output will be saved and choose file name

parser.add_argument("-a","--Regions", help="Write the path to the Bed file")
parser.add_argument("-b","--Reads", help="Write the path to the Bam file")
parser.add_argument("-c","--X_Cut_off_list", nargs='+', help="Write coverage depths thresholds, the first value will be used for coverage depth analytics and all three for coverage breadth")
parser.add_argument("-o", "--output_folder_name", nargs='+', help = "Write the name of the folder where the data is saved followed by start of file names" )
parser.add_argument("-e", "--executedBy", help="Enter name of person that run the coverage tool")

# Additional optional input settings:
# filtering of reads, mapQ and Phred score (Q score)
# Own command to samtools   
# validation list, contains regions with coverage depth values below 95% at first threshold
# combine rows, combine ROI from the same for example gene
# Strand Specific reads are used two generate two additonal coverage breadth tables one for the reversed and one for the forward reads
# Figures, pie chart, bar plot and a regions plot. Obs! The bar plot is only generated if combineRows is activated
# Hotspots, add positions of interest in the region figure. These will be marked by an arrow in the region figure
# Low regions, add regions known to be low in the full and short mean list as an extra column known = Yes or No. Yes if known to be low. Marked with a red color in the region position plot
# Exon and Transcript information in the mean coverage lists
# Detailed coverage saves the per base coverage depth list as a file

parser.add_argument("-p", "--phred_score_mapQ", nargs='+', help="Choose filter option for the reads. By adding 'all' the whole analysis will be computed with the filtered reads. If not specified the filtered reads will only be used as an additional row in the statistics table. Next the phred score threshold is specified followed by the mapQ threshold. example usage -p all 20 10 or -p 20 10")
parser.add_argument("-i", "--ownInput", help="Write the command to be sent to samtools as a string")
parser.add_argument("-v", "--validation", action="store_true", help="turn validation list option on")
parser.add_argument("-k", "--combineRows", action="store_true", help="turn per gene option on, combines rows in the bedfile from the same gene for example")
parser.add_argument("-s","--strandSpecific", action ="store_true", help="Create an additional statistics table with forward and reverse read coverage calculated seperatly")
parser.add_argument("-f", "--figures", action="store_true", help="Create figures") 
parser.add_argument("-t", "--hotspot", help="Add arrows to the region plot figure to indicate positions of interest")
parser.add_argument("-l", "--lowRegions", help="List of known to be low regions")
parser.add_argument("-n", "--ExonTranscript", action="store_true", help="adds exon and transcript information to the low coverage mean region list")
parser.add_argument("-d", "--detailedCoverage", action="store_true", help="Outputs a per position coverage depth list of the region bed file")

args = parser.parse_args()

# The current program version
ProgramVersion = "Program Version 1.0"


########################### Check that all mandatory parameter input is given

if args.Regions == None or args.Reads == None or len(args.X_Cut_off_list) < 3 or args.executedBy== None or len(args.output_folder_name) < 2:
	sys.exit("Not enough input: example -a myRegions.beg -b myReads.bam -c 10 20 30 -e UserName -o Folder Name and file name. Use help for more information")

############################# RUN THE PROGRAM ###################################################

print("Running coverage analysis ...")


# Create the an output folder for the analytic report
command = "mkdir " + str(args.output_folder_name[0])
s = subprocess.Popen(command, shell=True)
s.communicate()

####################### Create per base pair coverage list ###################################################################

# Open the Region bed file and convert into a list
Regions = []
with open(str(args.Regions), "r") as myfile:
	for line in myfile:
		element = line.strip('\n').split('\t')
		Regions.append([element[0], element[1],  element[2], element[3]])
myfile.close()



# To run the program with non-filtered bam file or your own samtool command as main option
if args.phred_score_mapQ == None or not(str(args.phred_score_mapQ[0]) == 'all'):

	# Set the data type variable to raw data, this will be shown in the statistics table
	dataType = ''

	# Calculate the coverage depth with samtools depth, eihter with own command from the user or by the default command
	if args.ownInput:
		command = args.ownInput + " -b " + str(args.Regions) + " " + str(args.Reads) + " > " + str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + "_coverage.tsv"
	else:
		command = "samtools depth -a -d 30000 -b " + str(args.Regions) + " " + str(args.Reads) + " > " + str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + "_coverage.tsv"
	s = subprocess.Popen(command, shell=True)
	s.communicate()

	# Open the filtered coverage file generated with samtools depth and convert into a list
	Reads_coverage = []
	with open(str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + "_coverage.tsv", "r") as myfile:
		for line in myfile:
			element = line.strip('\n').split('\t')
			Reads_coverage.append([element[0], element[1], element[2]])
	myfile.close()

	# Generate a detailed list with sublists of coverage values from the same region in the bed file. 
	import Detailed_generator
	detailed_list = Detailed_generator.detail_samtools(Regions, Reads_coverage)


####################### Filter Reads by phredscore and mapQ (Optional) ############################################################  
# Use samtools to filter the bam file with the given threshold value. This can either be used for the whole analysis or only as an additional field in the statistics table. 
# With the input -p all 20 10, the reads will be filtered by only containg bases over phred score 20 and reads above mapQ 10 and these reads are used for all calclations. 
# While -p 20 10 will filter the bam file and add to the statistics table with the non filtered bam file as mainly used. 

detailed_list_filter=[]
phred_score = False
mapQ = 0

# Check the type of filtering setting, save the phred score and mapQ as variables
if args.phred_score_mapQ:
	
	if str(args.phred_score_mapQ[0]) == 'all':

		print("Creates phred score and mapQ filtered coverage values for the analysis ...")
		phred_score=args.phred_score_mapQ[1]
		mapQ = args.phred_score_mapQ[2]
		dataType= 'Filtered'

	else:

		print("Creates phred score and mapQ filtered coverage values for the statistics table ...")
		phred_score = args.phred_score_mapQ[0]
		mapQ = args.phred_score_mapQ[1]
		dataType=''


	# Calculate the coverage depth with samtools depth
	if args.ownInput:
		command = args.ownInput + " -b " + str(args.Regions) + " " + str(args.Reads) + " > " + str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + "_filtered_coverage.tsv"
	else:
		command = "samtools depth -a -d 30000 -b " + str(args.Regions) + " -q " + str(phred_score) + " -Q " + str(mapQ) + " " + str(args.Reads) + " > " + str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + "_filtered_coverage.tsv"
	s = subprocess.Popen(command, shell=True)
	s.communicate()

	# Open the filtered coverage file generated with samtools depth and convert into a list
	Reads_filtered = []
	with open(str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + "_filtered_coverage.tsv", "r") as myfile:
		for line in myfile:
			element = line.strip('\n').split('\t')
			Reads_filtered.append([element[0], element[1], element[2]])
	myfile.close()


	# Generate a detailed list with sublists of coverage values from the same region in the bed file. 
	import Detailed_generator
	if str(args.phred_score_mapQ[0])=="all":
		detailed_list = Detailed_generator.detail_samtools(Regions, Reads_filtered)
	else:
		detailed_list_filter = Detailed_generator.detail_samtools(Regions, Reads_filtered)


############################# Filter bam file with reads to only contain + or - strand (OPTIONAL) #################################################

# These strand specific detailed coverage lists will be used to create statistic tables containing coverage breadth values at the three threshold values
if args.strandSpecific:
	print("Creates strand specific bam files for the additional statistics table ...")

	command_1 = "samtools view -F 0x10 -b " + str(args.Reads) + " > " + str(args.output_folder_name[0]) + "/outStrandPos.bam"
	command_2 = "samtools view -f 0x10 -b " + str(args.Reads) + " > " + str(args.output_folder_name[0]) + "/outStrandNeg.bam"
	s1 = subprocess.Popen(command_1, shell=True)
	s1.communicate()
	s2 = subprocess.Popen(command_2, shell=True)
	s2.communicate()

	if args.ownInput:
		depthCommand1 = args.ownInput + " -b " + str(args.Regions) + " " + str(args.output_folder_name[0]) + "/outStrandPos.bam > " + str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + "_cov_positive.tsv"
		depthCommand2 = args.ownInput + " -b " + str(args.Regions) + " " + str(args.output_folder_name[0]) + "/outStrandNeg.bam > " + str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + "_cov_negative.tsv"
	else:
		depthCommand1 = "samtools depth -a -d 30000 -b " + str(args.Regions) + " " + str(args.output_folder_name[0]) + "/outStrandPos.bam > " + str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + "_cov_positive.tsv"
		depthCommand2 = "samtools depth -a -d 30000 -b " + str(args.Regions) + " " + str(args.output_folder_name[0]) + "/outStrandNeg.bam > " + str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + "_cov_negative.tsv"
	
	s1 = subprocess.Popen(depthCommand1, shell=True)
	s1.communicate()
	s2 = subprocess.Popen(depthCommand2, shell=True)
	s2.communicate()

	# Open the forward coverage file generated with samtools depth and convert into a list
	Reads_positive = []
	with open(str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + "_cov_positive.tsv", "r") as myfile:
		for line in myfile:
			element = line.strip('\n').split('\t')
			Reads_positive.append([element[0], element[1], element[2]])
	myfile.close()

	# Open the reverse coverage file generated with samtools depth and convert into a list
	Reads_negative = []
	with open(str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + "_cov_negative.tsv", "r") as myfile:
		for line in myfile:
			element = line.strip('\n').split('\t')
			Reads_negative.append([element[0], element[1], element[2]])
	myfile.close()


	# Create per position coverage for the strand specific bam files 
	detailed_list_positive = Detailed_generator.detail_samtools(Regions, Reads_positive)

	detailed_list_negative = Detailed_generator.detail_samtools(Regions, Reads_negative)

###################### CombineRows (OPTIONAL) ################################################################################################## 
# Merge rows in the per position coverage depth lists (detailed lists) so that ALL coverage values from 
# the same, for example, gene will be in the same row. This is essential for per gene calculations. The output will be a new formated 
# detailed coverage result list and a combined region name list. The region name in the bedfile must be seperated with a '.' .Example: gene1.exon.2 etc.

det_list_filtered_formated =[]
splice=[] # Gives the indices of the subregions within the combined rows, this is used for the region figure
Info=[] # Gives the chromosome, start and stop for each subregion in the combined row region

# Create a list that only containing the region name column
Regions_list=[]
for line in Regions:
	Regions_list.append(line[3])

# Combine rows in the bedfile that has the same region name before the first dot as separator. Example: the rows Gene1.Exon.2 and Gene1.Exon.3 will be merged as Gene1
if args.combineRows:
	import CombineRows_generator
	print("Combine rows ...")

	# Formats the detailed coverage list by combining rows
	RegionNames, detailed_list_formated, splice = CombineRows_generator.CombineRowsList(detailed_list, Regions_list)

	# If the filtering option is turned on the filtered detailed coverage list is fomrated by merging rows  
	if args.phred_score_mapQ and not(args.phred_score_mapQ[0]=='all'):
		PhredRegionNames, det_list_filtered_formated, splice_phred = CombineRows_generator.CombineRowsList(detailed_list_filter, Regions_list)

	if args.strandSpecific:
		# If strand specific option turned on the negative and positive strand coverage detailed lists are formated by merging rows
		PosRegionNames, detailed_list_positive_formated, splice_pos = CombineRows_generator.CombineRowsList(detailed_list_positive, Regions_list)
		NegRegionNames, detailed_list_negative_formated, splice_neg = CombineRows_generator.CombineRowsList(detailed_list_negative, Regions_list)



###################### Create mean coverage regions list ##############################################################################################
# In the Mean coverage list subregions including values all above the coverage threshold or all bellow is saved from the per base coverage list.
# From the mean coverage list only the low covarge subregions is saved to a new short mean coverage list. 
# Example: detailed coverage list = [1, 1, 1, 4, 5, 2, 2] and if the coverage threshold is = 3 
# 	=> Mean coverage list = [[start = 0, stop = 3, mean = 1], [3, 5, 4.5], [5, 7, 2]]
# 	=> The short mean list = [[0,3,1],[5,7,2]] 

print("Generating Mean coverage lists ...")


# create a list with only start positions from each line in the bedfile,
# this list is used to compute new start and stop positions for the new subregions
Start_pos_list =[]
for line in Regions:
	Start_pos_list.append(int(line[1]))

## Add columns to the mean region list
info_list =[] 
info_temp=[]
mean_index=0

if args.ExonTranscript: # Adds exon number, transcript and chromosome. Extract information from the region name
	for element in Regions_list:
		info_temp = element.split('.')
		info_temp = [info_temp[2], info_temp[4], info_temp[5]]
		info_list.append(info_temp)
		mean_index=5

else: # Only add the chromosome not exon and transcript, extract chromosome from the bed file column
	for element in Regions:
		info_temp = element[0]
		info_list.append([info_temp])
		mean_index=3


# Create the Mean coverage list
import Meanlist_generator
Row_temp = []
MeanCov_List = []
index = 0

################ Create the full coverage mean list
# Compute the mean sub regions for each row in the detailed list and create the Mean Coverage list 

for line in detailed_list:
	Row_temp = Meanlist_generator.Meanfunction(line, args.X_Cut_off_list[0], Start_pos_list[index], info_list[index])
	MeanCov_List.append(Row_temp)
	index += 1

# Merge the region rows in the mean coverage list generated above, so that regions from the same for example gene is merged. Only if combinerows is activated
if args.combineRows:
	import CombineRows_generator
	MeanNames, MeanCoverage_formated, splice_meanList = CombineRows_generator.CombineRowsList(MeanCov_List, Regions_list) 


################# Create a short list of mean coverage regions under the coverage depth threshold

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


################################# Known low regions as input (OPTIONAL) #########################################

if args.lowRegions:
	#Extract the name and information of the low regions
	name_lowRegions=[]
	lowRegions=[]
	with open(str(args.lowRegions), "r") as myfile:
		for line in myfile:
			element = line.strip('\n').split('\t')
			name_lowRegions.append(element[3])
			lowRegions.append([element[0], element[1], element[2], element[3]])

	myfile.close()

	# Save the indices of the low regions this will be added to the mean short list and the region figure
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
					if str(line[mean_index-3]) == str(Low_region_compare[0]) and int(line[mean_index-2])==int(Low_region_compare[1]) and int(line[mean_index-1])==int(Low_region_compare[2]):

						LowRegionIndex.append([index, count])
					count+=1
				count=0	
		index+=1

	# Add the low regions as a column in Mean coverage short list

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
	hotspots=[]

	# Extract the name and information of the hotspots

	with open(str(args.hotspot), "r") as myfile:
		for line in myfile:
			element = line.strip('\n').split('\t')
			name_hotspots.append(element[3])
			hotspots.append([element[0], element[1], element[2], element[3]])
	myfile.close()


	#Calculate the new start positions of the hotspots that will be used in the region figure. The new start indices is the number of bases before the hotspot in the same region
	if args.combineRows:
		start_positionlist=[]
		length_sum=0
		index=0
		for hline in hotspots:
			for region in Regions:
				if str(region[3].split('.')[0])==str(hline[3]): 

					if int(hline[2]) > int(region[2]):
						length_sum+=(int(region[2])-int(region[1]))
					if int(hline[1]) > int(region[1]) and int(hline[1]) <= int(region[2]):
						length_sum+=int(hline[2])-int(region[1])

			start_positionlist.append(length_sum)
			length_sum=0
			index+=1

	else:
		start_positionlist=[]
		length_sum=0
		index=0
		for hline in hotspots:
			for region in Regions:
				if str(region[3])==str(hline[3]): 

					if int(hline[2]) > int(region[2]):
						length_sum+=(int(region[2])-int(region[1]))
					if int(hline[1]) > int(region[1]) and int(hline[1]) <= int(region[2]):
						length_sum+=int(hline[2])-int(region[1])

			start_positionlist.append(length_sum)
			length_sum=0
			index+=1


	# Save the indices of the hotspots this will be used for the region figure, [row index, index of subregion, start index in region]  
	start_index=0
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
					if str(line[mean_index-3]) == str(hotspots_compare[0]) and int(line[mean_index-2]) <= int(hotspots_compare[1]) and int(line[mean_index-1]) >= int(hotspots_compare[2]):

						hotspotIndex.append([index, count, start_positionlist[start_index]])
						start_index+=1

					count+=1

				count=0
		index+=1

			 

###################### Create statistics table #####################################################################################################
# Uses the detailed coverage list and reports coverage breadth at and over the choosen coverage depth thresholds. 
# The statistics table can also contain additional results from a filtered bam file and an optional validation column. 
# The validation option checks if 95% coverage bredth or more at the first coverage threshold value. If under 95% the column is marked with **** and the
# region is added to the validation list. If the strand specific option is turned on an additional statistics tables with strand specific reads are generated.

print("Generating Statistics table ...")

# If validation option is turned on, create a list with region info for the validation table. With the columns chr, start, stop, length
Region_info=[]

if args.validation:
	for line in Regions:
		Region_info.append([line[0], line[1], line[2], int(line[2])-int(line[1])])

	# If combine rows activated, remove duplicates in the regon info list
	if args.combineRows:
		import CombRegionInfo
		Region_info = CombRegionInfo.CombineRegionInfo(Region_info, Regions_list)

# Decide on what coverage lists to use for the statistics table
if args.combineRows:
	detailed_list_stat = detailed_list_formated
	detailed_filt_stat = det_list_filtered_formated

	if args.strandSpecific:
		detailed_positive = detailed_list_positive_formated
		detailed_negative = detailed_list_negative_formated
	Stat_table_names = RegionNames
else:
	detailed_list_stat = detailed_list
	detailed_filt_stat = detailed_list_filter

	if args.strandSpecific:
		detailed_positive = detailed_list_positive
		detailed_negative = detailed_list_negative
	Stat_table_names = Regions_list

	
import Statistics_generator

# Compute statistics table 	
stat_table, stat_table_filter, validation_list = Statistics_generator.stat_table(detailed_list_stat, Stat_table_names, args.validation, phred_score, detailed_filt_stat, args.X_Cut_off_list, Region_info, dataType)    

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
Logg.append(["Region file (BED file): ", Region_fileName])
Logg.append(["Read file (BAM file): ", Reads_fileName])
Logg.append(["Output folder and file name:", args.output_folder_name[0], args.output_folder_name[1]])
Logg.append([" "])

# Add the chosen coverage depth thresholds
Logg.append(["Coverage depth thresholds: ", args.X_Cut_off_list[0], args.X_Cut_off_list[1], args.X_Cut_off_list[2]])

# Calculate the mean coverage breadth for the columns in the statistics table and add to the logg file
stat_temp1 = 0.0
stat_temp2 = 0.0
stat_temp3 = 0.0

for line in stat_table:
	stat_temp1 += line[2]
	stat_temp2 += line[3]
	stat_temp3 += line[4]

Logg.append(["Mean Coverage Breadth: ", round(float(stat_temp1)/float(len(stat_table)),2),round(float(stat_temp2)/float(len(stat_table)),2), round(float(stat_temp3)/float(len(stat_table)),2)])

#Calulate the total mean coverage and add to the log file
Total_Mean_Cov = 0
count=0

for line in Mean_full_list: 
	for element in line:
		Total_Mean_Cov += element[mean_index]
		count+=1

Total_Mean_Cov = round(float(Total_Mean_Cov)/float(count),2)
Logg.append(["Mean Coverage Depth: ", str(Total_Mean_Cov) + ' X'])

# Optional input:

Logg.append([" "])

# filter option of bam file
if args.phred_score_mapQ:
	if str(args.phred_score_mapQ[0])=="all":
		Logg.append(["Phred score: ", phred_score, "mapQ: ", mapQ, "Used in all calculations"])
	else:
		Logg.append(["Phred score: ", phred_score, "mapQ: ", mapQ, "Used only in an additional row in the statistics table"])

if args.ownInput:
	Logg.append(["Samtools command: ", args.ownInput])

# validation option
if args.validation:
	Logg.append(["Validation list generated"])

# Combine rows
if args.combineRows:
	Logg.append(["Combine regions activated"])

# Generation of figures
if args.figures:
	Logg.append(["Figures generated"])
	# If hotspots been added to the region figure
	if args.hotspot: 
		Logg.append(["Hotspots added"]) # add what hotspot list that was used

# If exon and transcript has been added to the mean coverage list
if args.ExonTranscript:
	Logg.append(["Exon and transcript added"])

# If a list of low regions was given and an additional column was added to the mean coverage list
if args.lowRegions:
	Logg.append(["Low region list added", args.lowRegions])

# If the detailed list of coverage depth per position was saved as output
if args.detailedCoverage:
	Logg.append(["Per base coverage depth list generated"])

# If strand specific reads are generated
if args.strandSpecific:
	Logg.append(["Strand specific tables generated"])

	


######################################## Figures ######################################################################################################
if args.figures == True:

	# Choose what detailed coverage list and region name list to use
	if args.combineRows:
		list_for_figures = detailed_list_formated
		Names_for_figures = RegionNames

	else:
		list_for_figures = detailed_list
		Names_for_figures=Regions_list
		Names_for_regionfig=[]

		for name in Regions_list:
			Names_for_regionfig.append(name.split('.')[0])


	##### PIE CHART #######################################################################################################################################
	print("Generating figures ...")

	# One pie chart per row in bedfile (region file) or one per combined region. The pie charts are only created for regions that has some subpart below the coverage depth threshold, X_cut_off_list[0]
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	import tkinter
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

		# Only create pie charts with lines that has some position under the coverage threshold
		if not(int(cov_under) == 0):	
			fig = PieChart.Generate_Pie(line, args.X_Cut_off_list[0], index_order, fig, Names_for_figures[Name_index], text_size)

			# when 4 pie charts exists in the figure they are saved in a pdf file
			if index_order==4:
				# Save the pie charts as a PDF 
				fig.savefig(str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + '_Cov_pie_' + str(image_count) + '.pdf')
				fig = plt.figure()
				index_order=0
				image_count+=1
			index_order+=1
		cov_under=0
		Name_index+=1

	# Make sure that the last image is printed even if there is less than 4 pie charts in it
	if not(fig == None): 
		fig.savefig(str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + '_Cov_pie_' + str(image_count) + '.pdf')

############################################# Bar plot

	if args.combineRows:

		# Create a list with subpart names that is used as x-axis labels for each bar
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
					fig.savefig(str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + '_Cov_bar_' + str(image_count) + '.pdf')
					fig = plt.figure() # Create a new figure to draw the pie charts in
					index_order=0
					image_count+=1
				index_order+=1
			Name_index+=1
			bellow=0

		# Make sure that the last image is printed even if there is less that 4 pie charts in it
		if not(fig == None): 
			fig.savefig(str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + '_Cov_bar_' + str(image_count) + '.pdf')

	###### Per region coverage plot ###############################################################################################################


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
	only_under_index = 0
	data = []
	column =[]
	rows=[]
	fig = plt.figure()
	splice_index=0
	hotspots_arrows=[]

	if splice == None: 
		splice=[]

	# Create a region figure for all regions that has coverage values below the threshold
	for line in list_for_figures:

		for element in line: 

			if int(element) < int(args.X_Cut_off_list[0]):
				count_under+=1

		if int(count_under) > 0:

			# if the hotspot option is active find all hotspots for the current region and save the start position of those in hotspots_arrow
			if args.hotspot:
				for hs in hotspotIndex:
					if int(hs[0]) == int(Name_index):
						hotspots_arrows.append(int(hs[2]))

			# Generate the figure, with splice if combined rows otherwise splice is set to false
			if args.combineRows:
				fig = Region_position_plot.Region_generator_plot(line, Names_for_figures[Name_index], args.X_Cut_off_list[0], index_order, fig, splice[splice_index], hotspots_arrows)
			else:
				fig = Region_position_plot.Region_generator_plot(line, Names_for_regionfig[Name_index], args.X_Cut_off_list[0], index_order, fig, False, hotspots_arrows)
			hotspots_arrows=[]
			index_order+=1

			# The images is saved as a PDF with 2 regions figures and 2 tables next to the figures, 4*4 image
			ax = fig.add_subplot(2,2,index_order)

			# Append the data for the region figure table
			for element in Text_list[only_under_index]:
				data.append(element)

			# Define the labels in the table
			if not(args.ExonTranscript):	
				column = ['Chr','Start', 'Stop', 'Mean','Length']
			else:
				column = ['Exon', 'Transcript', 'Chr', 'Start', 'Stop', 'Mean', 'Length']

			# Name the low regions, first region is called R1 next R2 etc. 
			for i in range(len(data)):
				rows.append('R' +  str(i+1))

			# If low regions are sent as input, add red color to any line in the table containg a known to be low region
			color=[]
			lowregions_conter=0

			if args.lowRegions:
				for row_data in data:
					for low_row in lowRegions:
						if str(row_data[mean_index-3]) == str(low_row[0]) and int(row_data[mean_index-2]) == int(low_row[1]) and int(row_data[mean_index-1]) == int(low_row[2]):
							lowregions_conter=1
					if lowregions_conter==1:
						lowregions_conter==0
						color.append('tomato')
					else:
						color.append('white')

			if color==[]:
				fig_table = ax.table(cellText= data, rowLabels = rows, colLabels = column, loc='center')
			else:
				fig_table = ax.table(cellText= data, rowLabels = rows, colLabels = column, loc='center', rowColours=color)
			data=[]
			rows=[]
			color=[]
			ax.axis('off')

			# Save the figures to a PDF		
			if index_order == 4:
				fig.savefig(str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + 'RegionCoverage_' + str(image_count) + '.pdf')
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
		fig.savefig(str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + 'RegionCoverage_' + str(image_count) + '.pdf')

	

################## Save list and tables to csv files ###################################################################################################

print("Saving to csv files ...")

###### Save the non strand specific statistics table to a csv file

index=0
with open(str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + '_Stat_table.csv', 'w') as myfile:
	wr = csv.writer(myfile)

	if args.validation:
		wr.writerow(['Data type','Region Name', str(args.X_Cut_off_list[0]) +'X', str(args.X_Cut_off_list[1])+'X', str(args.X_Cut_off_list[2])+'X', 'Validation'])
	else:
		wr.writerow(['Data type','Region Name', str(args.X_Cut_off_list[0]) +'X', str(args.X_Cut_off_list[1])+'X', str(args.X_Cut_off_list[2])+'X'])
	for row in stat_table:
		wr.writerow(row)
		if args.phred_score_mapQ and not(str(args.phred_score_mapQ[0])) == "all":
			wr.writerow(stat_table_filter[index])
			index+=1
myfile.close()

###### Save the positive strand specific statistics table to a csv file
if args.strandSpecific:
	index=0
	with open(str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + '_Stat_table_positiveStrand.csv', 'w') as myfile:
		wr = csv.writer(myfile)
		wr.writerow(['Region Name', str(args.X_Cut_off_list[0]) +'X', str(args.X_Cut_off_list[1])+'X', str(args.X_Cut_off_list[2])+'X', 'Strand difference'])
		for row in Stat_table_positive_final:
			wr.writerow(row)
	myfile.close()

###### Save the negative strand specific statistics table to a csv file

	index=0
	with open(str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + '_Stat_table_negativeStrand.csv', 'w') as myfile:
		wr = csv.writer(myfile)
		wr.writerow(['Region Name', str(args.X_Cut_off_list[0]) +'X', str(args.X_Cut_off_list[1])+'X', str(args.X_Cut_off_list[2])+'X', 'Strand difference'])
		for row in Stat_table_negative_final:
			wr.writerow(row)
	myfile.close()



###### save Mean coverage short list to csv file

index=0
with open(str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + '_MeanCoverageShortList.csv', 'w') as myfile:
	wr =csv.writer(myfile)
	if args.ExonTranscript:
		if args.lowRegions:
			wr.writerow(['Region Name', 'Exon', 'Transcript', 'Chr', 'Start', 'Stop', 'Mean Coverage', 'Length', 'Known'])
		else:
			wr.writerow(['Region Name', 'Exon', 'Transcript', 'Chr','Start', 'Stop', 'Mean Coverage', 'Length'])
	else:
		if args.lowRegions:
			wr.writerow(['Region Name', 'Chr', 'Start', 'Stop', 'Mean Coverage', 'Length', 'Known'])
		else: 
			wr.writerow(['Region Name', 'Chr', 'Start', 'Stop', 'Mean Coverage', 'Length'])

	for row in MeanCov_shortList:
		if not(row == []):
			if args.lowRegions:
				line = [MeanNames[index]] + Mean_Short_with_lowRegionInfo[index]
				wr.writerow(line)
			else:
				wr.writerow([MeanNames[index]] +row)

		index+=1
myfile.close()

###### save Mean Coverage Full List to csv file

index=0
rowtemp=[]
index_2=0

with open(str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + '_MeanCoverageFullList.csv', 'w') as myfile:
	wr= csv.writer(myfile)

	wr.writerow(['Region Name', 'Chr', 'Start', 'Stop', 'Mean Coverage', 'Length'])

	for row in MeanCov_List:

		if int(len(row)) > 1:

			for element in row:
				rowtemp+=row[index_2]
				index_2+=1
			wr.writerow([Regions_list[index]] + rowtemp)
			rowtemp=[]
			index_2=0

		else:	
			wr.writerow([Regions_list[index]] + row[0])
			
		index+=1


myfile.close()

##### Save the log file to csv file

with open(str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + '_Log.csv', 'w') as myfile:
	wr = csv.writer(myfile)
	wr.writerow(['LOG'])
	for line in Logg:
		wr.writerow(line)
myfile.close()

# Save the validation list for the original statistics table to a csv file

if args.validation and not(validation_list == None):

	with open(str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + '_Validation_list.csv', 'w') as myfile:
		wr = csv.writer(myfile)
		wr.writerow(['Region Name', '% Coverage at: ' + str(args.X_Cut_off_list[0])+ ' X', 'Chr', 'Start', 'Stop', 'Length']) 
		for line in validation_list:
			wr.writerow(line)
	myfile.close()

# Saves the per base coverage depth list to a csv file

if args.detailedCoverage:
	index=0
	Region_info_det=[]
	for line in Regions:
		Region_info_det.append([line[3], line[0], line[1], line[2]])

	with open(str(args.output_folder_name[0]) + "/"+ str(args.output_folder_name[1]) + '_PerBaseCoverage.csv', 'w') as myfile:
		wr =csv.writer(myfile)
		wr.writerow(['Region Name', 'Chr', 'Start', 'Stop','Coverage Depth'])
		for line in detailed_list: 
			line = Region_info_det[index] + line
			wr.writerow(line)
			index+=1
	myfile.close()
