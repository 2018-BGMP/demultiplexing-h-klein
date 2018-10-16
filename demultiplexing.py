# Helena Klein
#Demultiplexing Part 2
# To use: 
# local use: python demultiplexing.py -R1 Unit_Test/read_test_1.txt -R2 Unit_Test/read_test_2.txt -I1 Unit_Test/index_test_1.txt -I2 Unit_Test/index_test_2.txt -key indexes_key.txt -T 30
#Talapas Usage: python3 demultiplexing.py -R1 <read1.fq> -R2 <read2.fq> -I1 <index1.fq> -I2 <index2.fq> -T <mean read threshold>
# RUN IN OUTPUT DIRECTORY

# This approach only works if the reads in all of the files are in the same order. But if they are in the same order
# Across the read and index files, this should avoid iterating multiple times through any of the files and saving
# a ridiculous amount to memory. 

import argparse
import gzip

def get_arguments():
	parser = argparse.ArgumentParser(description="Demultiplex files and report index hopping given paired end read files, their index files, and known indexes")
	parser.add_argument('-R1', '--read_1', required=True, help='Indicate your first read file to demultiplex')
	parser.add_argument('-R2', '--read_2', required=True, help='Indicate your second read file to demultiplex')
	parser.add_argument('-I1', '--index_1', required=True, help='Indicate your first index file for use in demultiplexing and measureing index hopping')
	parser.add_argument('-I2', '--index_2', required=True, help='Indicate your second index file for use in demultiplexing and measureing index hopping')
	parser.add_argument('-T', '--threshold', required=True, help='Indicate your mean read quality threshold for quality filtering')
	parser.add_argument('-key', '--index_key', required=True, help='Provide a txt file with the expected indexes')

	return parser.parse_args()
	
args = get_arguments()
R1 = args.read_1
R2 = args.read_2
I1 = args.index_1
I2 = args.index_2
KEY = args.index_key
threshold = float(args.threshold)

def convert_phred(letter):
    """Converts a single character into a phred score"""
    phred = ord(letter)-33
    
    return phred

# Open all generated index files using a funciton

def open_output(dictionary, read_number):
	"""Opens all files that are keys in a dictionary when given a dictionary, and generates a dictionary of file handles for future use"""
	handle_dict = {}
	for index in dictionary:
		handle_format = 'R{}_{}'
		handle = handle_format.format(read_number, index)
		handle_dict[index] = handle
	for index in dictionary:
		handle_dict[index]= open(dictionary[index], 'a')
        
    
	return handle_dict

# close all generated index files using a function 
def close_output(dictionary):
	"""Closes all files in the output dictionary from 'open_output'"""
	for index in dictionary:
		dictionary[index].close()
    
	return    

# This function takes the reverse complement of whatever I stick in 
def reverse_complement(sequence):
	"""Gives the reverse complement of the input sequence"""
	complement = {'A':"T", 'G':'C', 'C':'G', 'T':'A'}
	reverse = ""
	for letter in sequence:
		reverse = complement[letter]+reverse

	return reverse
    
# Make all the output file names for each index using a file with all of the indexes listed
# This section also generates a dictionary to count the appearance of each properly matched
# index

with open(KEY, 'rt') as key:
	F_1 = 'R1_{}.fq'
	F_2 = 'R2_{}.fq'
	R1_output={}
	R2_output={}
	index_counts = {}
    
	for line in key:
		line=line.strip()
		R1_output[line] = F_1.format(line)
		R2_output[line] = F_2.format(line)
		index_counts[line]=0
        
# Put saved items in this format when writing to the different read formats.
read_to_write='{}:{}\n{}\n+\n{}\n'


# Open both index files at once and compare across them, and use that to open both the read files at once
with gzip.open(I1, 'rt') as i1, gzip.open(I2, 'rt') as i2, \
gzip.open(R1, 'rt') as r1, gzip.open(R2, 'rt') as r2:    
	LN = 0
	undetermined = 0
	swapped = 0
	matched = 0
	low_quality = 0
	total_reads=0
	output_counts = {}
        
	for line1, line2, read1, read2 in zip(i1, i2, r1, r2):
		LN += 1
		read1=read1.strip()
		read2=read2.strip()
		line1 = line1.strip()
		line2 = line2.strip()
		if LN%4==1:
			#counts header lines for final read count
			total_reads += 1
			R1_head = read1
			R2_head = read2
                
# This section compares the indexes directly as we iterate through them
		elif LN%4==2:
		#Holds on to the read sequence lines for later
			R1_sequence = read1
			R2_sequence = read2
			index = line1
			# save the index for appending to header later
			index_R1 = line1
			index_R2 = line2

			#or if one contains an N    
			if "N" in line1+line2:
				result='undetermined'
				undetermined += 1
			elif line1==reverse_complement(line2):
				result = 'match'
				if line1 not in output_counts:
					output_counts[line1] = 1
				else:
					output_counts[line1] +=1
			else:
				result="swapped"
				swapped += 1
                    
# This section should take the results from comparing indexes and put the full read into the correct output files    
		elif LN%4==0:
			R1_quality = read1
			R2_quality = read2
			if result == 'match':
			# This allows us to set a quality threshold for the index reads even if they are matched so that
			# we can filter out the low quality indexes. 
				read_sum = 0
				for position in index_R1+index_R2:
					index_sum += convert_phred(position)
				mean = index_sum/202
                    
				if mean > threshold:
					# This is the last check to make sure the indexes are matched
					matched += 1
                    
                    # Open files to write sorted outputs to
					R1_files= open_output(R1_output,1)
					if index in R1_files:
						R1_files[index].write(read_to_write.format(R1_head, index_R1, R1_sequence, R1_quality))
						
					# This removes low quality indexes that do not match the provided ones
					else:
						low_quality += 1
						R1_undetermined = open('R1_undetermined.fq', 'a')
						R1_undetermined.write(read_to_write.format(R1_head, index_R1, R1_sequence, R1_quality))
                    
					R2_files=open_output(R2_output, 2)
					if index in R2_files:
						R2_files[index].write(read_to_write.format(R2_head, index_R2, R2_sequence, R2_quality))
					# This removes low quality indexes that do not match the provided ones
					else:
						low_quality += 1
						R2_undetermined = open('R2_undetermined.fq', 'a')
						R2_undetermined.write(read_to_write.format(R2_head, index_R2, R2_sequence, R1_quality))
				else: 
					result = "low_quality"
					low_quality += 1
                        
					R1_undetermined = open('R1_undetermined.fq', 'a')
					R1_undetermined.write(read_to_write.format(R1_head, index_R1, R1_sequence, R1_quality))
					
					R2_undetermined = open('R2_undetermined.fq', 'a')
					R2_undetermined.write(read_to_write.format(R2_head, index_R2, R2_sequence, R1_quality))
                                   
			else:
                #write non-matches to "undetermined" file
				R1_undetermined = open('R1_undetermined.fq', 'a')
				R1_undetermined.write(read_to_write.format(R1_head, index_R1, R1_sequence, R1_quality))
                    
				R2_undetermined = open('R2_undetermined.fq', 'a')
				R2_undetermined.write(read_to_write.format(R2_head, index_R2, R2_sequence, R2_quality))
                #index.write(read_to_write.format(R1_head, R1_sequence, result, R1_quality))

close_output(R1_files)
close_output(R2_files)
R1_undetermined.close()
R2_undetermined.close()

swapping_output = 'Percentage of index swapping: {}%'


print(swapping_output.format(round((swapped/total_reads*100),2)))
print("Total Undetermined Reads:", undetermined)
print("Total Low-quality Reads:", low_quality)
print("Total Matched Reads:", matched)
print("Total Swapped Reads:", swapped)
print("Total Reads:", total_reads)
print()
print("Percentage of reads per index")
print("Index\tPercent of Total Sequencing Run")
for index in output_counts:
    output_counts[index]=round((output_counts[index]/total_reads)*100, 2)
    print(index, output_counts[index], sep='\t')

