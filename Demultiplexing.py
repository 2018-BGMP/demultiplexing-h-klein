# Helena Klein
#Demultiplexing Part 1
# To use: run in folder where I want the output images

import argparse
import numpy as np
#import matplotlib
#import matplotlib.pyplot as plt
import gzip

def get_arguments():
	parser = argparse.ArgumentParser(description="Plot the average quality score per base position in reads or indexes")
	parser.add_argument('-f', '--file', help='Indicate a file to average the quality scores')
	return parser.parse_args()
	
args = get_arguments()
file = args.file

files = {'/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz': 'read', '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz':'index', '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz':'read', '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz':'index', 'index_test.fq.gz': 'index', 'read_test.fq.gz': 'read'}


def convert_phred(letter):
    """Converts a single character into a phred score"""
    phred = ord(letter)-33
    
    return phred
    

def populate_array(file):
    '''This function loops through a FASTQ file and converts the Phred quality score from a letter to its corresponding number. This function adds it to an ongoing sum of quality scores for each position'''
    LN = 0
    number_quality = 0
    if files[file] == 'read':
        scores = np.zeros(shape=(101,1))
    elif files[file] == 'index':
        scores = np.zeros(shape=(8,1))
    
    with gzip.open(file, 'rt') as fh:
        for line in fh:
            line = line.strip()
            LN +=1
            
            if LN%4 == 0:
                number_quality += 1
                for position in range(len(line)):
                    scores[position] += convert_phred(line[position])
        
    return scores, LN, number_quality


# Run populate_array function in this cell
scores_sum, NR, NQ = populate_array(file)

# Calculating the mean quality score 
mean_scores = scores_sum/NQ

#Plotting!

y = mean_scores.tolist()
x = range(len(mean_scores))

#plt.figure(figsize=(16, 10), dpi= 80, facecolor='w', edgecolor='k')

#plt.plot(x, y, marker='o')
#plt.xlabel('Base Position')
#plt.ylabel('Mean Quality Score')
#plt.title('Mean Quality Scores from Illumina Sequencing')

#plt.show()

name = 'mean_{}.tsv'
filename = file.strip().split('/')[-1]
out = name.format(filename)

with open(out, 'w') as of:
	# Column 1 is x, column 2 is y
	for item in range(len(y)):
		line = str(item) + '\t' + str(y[item][0])+ '\n'
		of.write(line)
		
#plt.savefig(image_name.format(file))



