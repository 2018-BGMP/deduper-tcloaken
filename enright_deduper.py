#!/usr/bin/env python3
"""
author: Trevor Enright
date: Fall 2018
input SAM file

returns a sam file without PCR duplicates

assume genome reference is forward stranded "to the right"
forward strand is the same direction so we need to "left align"
reverse strand is opposite thus we need to "right align"
								   
								   
REF     :-------------------------->	
FOR_strand       :-->
REV_strand             <--:

umi list file:
umi1\n
umi2\n
umi3\n
...

for getting a list of UMIs in the command line:

cat <SAMfile> | grep -v ^@ | cut -f 1 | cut -d ":" -f<UMI field> | sort | uniq -c | sort -nr | head -n <top umis> | sed -E "s/\s\s[0-9]+\s//" > umi_list
"""


####################
# Import libraries #
####################
import sys
import argparse as ap
import re
import gzip

#import numpy as np


#####################
# Argument Parser   #
#####################
def get_args():
	parser = ap.ArgumentParser(description="This is a deduper for sorted sam files. (see samtools sort). This algorithm takes into account soft clipping, strandedness (requires strand specific data), and Unique molecular identifiers (UMIS), also deals with indels and splicing in CIGAR string. Does not currently support paire-end.  Does support UMI file or randomers.  Use [--strict n] option to control hamming_distance correction for umis/randomers")
	parser.add_argument("-f","--file", help="required arg, put the path the SAM file after the flag"
	,type=str,required=True)
	parser.add_argument("-p",dest='paired_end', action='store_true',help="optional arg, designates file is paired end (not single-end)"
	)
	parser.add_argument("-u","--umi",help="optional arg, designates file containing the list of UMIs (unset if randomers instead of UMIs)"
	,type=str)
	parser.add_argument("-s","--strict",help="optional arg, designates the restriction for hamming distance of UMIs to one another, accepts less than or equal to a number, default is 1"
	,type=int,default=0)
	
	return(parser.parse_args())

######################
# GLOBAL VARIABLES 1 #
######################
args = get_args()
file = args.file	#SAM file
UMIS = args.umi    #file with list of umis
paired = args.paired_end 
strict = args.strict

#if(paired)
#	sys.exit("Does not accept paired-end at this time") 
umi_list = []  #get the list of umis
if UMIS is not None:
	with open(UMIS,"r") as fh:
		for line in fh:
			umi_list.append(line.strip())
UMI_list_flag = (umi_list != []) #umi list is true if not empty


################
# FUNCTIONS    #
################

def hamming_distance(s1, s2):
    #Return the Hamming distance between equal-length sequences
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def CIGAR_FOR_PARSE(cigar_string):
	##input cigar string
	#returns number (int) to be added to position to get "left"-based
	#alignment
	match = re.search(r"^(\d+)S",cigar_string)
	if match is not None:
		return int(match.group(1)) #return the number \d
	else:
		return 0
		

def CIGAR_REV_PARSE(cigar_string):
	#input cigar string
	#returns number (int) to be added to position to get "right"-based
	#alignment, Adds right side soft clipping, Matches, Deletions, and N's
	cigar_string += "\n" #add end of line character
	match = re.search(r"(\d+)S\n",cigar_string)
	S = 0
	if match is not None:
		S = int(match.group(1)) #gets last soft clipping
	sum_nums = sum([int(i) for i in re.findall(r"(\d+)[DMN]",cigar_string)]) #gets digits for deletions, matches, and gaps
	
	X = sum_nums + S
	return X

	
def add_read(dict,read,positionY,umi,chrm):
	#input dictionary, read line, and reference Y position
	dict[(umi,chrm,positionY)] = read
	return None

def tidy_FOR_dict(dict,current_position,read_length):
	#We'll use the following logic to remove
	#extraneous reads in the Forward strand dictionary
	#"left aligned" positions < current_position - read_length
	# returning a new dictionary without the above positions
	
	new = {k:v for k,v in dict.items() if k[2] >= (int(current_position) - read_length)}		
	return new
	
def tidy_REV_dict(dict,current_position,read_length):
	#We'll use the following logic to remove
	#extraneous reads in the Reverse compliment strand dictionary
	#"right aligned" positions < current_position + read_length
	# returning a new dictionary without the above positions
	
	new = {k:v for k,v in dict.items() if k[2] >= (int(current_position) + read_length)}		
	return new	
		
def is_Duplicate(dict,umi,chrm,Y):
	#given a dictionary to search through for the umi,chrom,and position of interest
	#return a true or false,  True if in dictionary
	#first check if Y reference is in dictionary
	if (umi,chrm,Y) not in dict:
		#check if umi hamming_distance 
		#then check umi hamming_distance with all others that match chrm and Y
		if not UMI_list_flag:
			if strict > 0:
				ham_list = sorted([hamming_distance(x[0],umi) for x in dict if x[1]==chrm and x[2] == Y])
				if len(ham_list) == 0:
					#not a duplicate
					return False
				elif (len(ham_list) == 1):
					return (ham_list[0] <= strict)
				elif (ham_list[0] != ham_list[1]):
					#There's theres a minimum that matches and 
					#there's not anoter of the same hamming distance
					#we'll call this a duplicate!
					return (ham_list[0] <= strict)
				else:
					return False
			else:
				return False
		else:
			#not a duplicate
			return False
	else:
		#is a duplicate
		return True

def UmiListCorrection(UMI):
	"""
	only should be called if a list of UMIs is provided
	returns 3 tuple: a string (corrected umi or true umi)
	               and a bool (true if meets critrea/false if otherwise)
				   and a bool (true if original UMI not in UMI list,  False if otherise)
	"""
	
	if (UMI not in umi_list):
			#try to repair UMI?
			if strict > 0:					
				try:
					UMI_correction = umi_list[sorted([(hamming_distance(x,UMI),i) for i,x in enumerate(umi_list) if hamming_distance(x,UMI) <= strict])[0][1]]
					return (UMI_correction,True,True)
				except:
					return (UMI,False,True) #couldn't find correction
			else:
				#not in list and no correction (skip)
				return (UMI,False,True)
	else:
		return (UMI,True,False)
	

def add_paired_read(dict,read1,read2,positionReverse,positionForward,umi1,umi2,chrm,flag):
	#input dictionary, read lines positionally
	if ((int(read1.split("\t")[1]) & 64) == 64): #read one first?
		dict[(umi1,umi2,chrm,positionForward,positionReverse,flag)] = read1+read2
	else:
		dict[(umi1,umi2,chrm,positionForward,positionReverse,flag)] = read2+read1
	return None			
	
def isPairedDuplicate(dict,umi1,umi2,chrm,forwardA,reverseA,flag):
	#returns if PCR duplicate or not (bool True if duplicate, False if otherwise)
	if (umi1,umi2,chrm,forwardA,reverseA,flag) not in dict:
		#check if umi hamming_distance 
		#then check umi hamming_distance with all others that match chrm and Y
		if not UMI_list_flag: #if there's not a list of umis
			if strict > 0:
				ham_list1 = sorted([hamming_distance(x[0],umi1) for x in dict if x[2]==chrm and x[3] == forwardA and x[4] == reverseA] and x[5] == flag)
				ham_list2 = sorted([hamming_distance(x[1],umi2) for x in dict if x[2]==chrm and x[3] == forwardA and x[4] == reverseA] and x[5] == flag)
				if len(ham_list1) == 0 or len(ham_list2) == 0:
					#not a duplicate
					return False
				elif len(ham_list1) == 1 or len(ham_list2) == 1:
					return (ham_list1[0] <= strict and ham_list2[0] <= strict)
				elif ham_list1[0] != ham_list1[1] and ham_list2[0] != ham_list2[1]:
					return (ham_list1[0] <= strict and ham_list2[0] <= strict)
					
			else: # can't correct for UMI so its not a duplicate by default
				return False
		else:
			#not a duplicate, since theres a list of umis and they are already adjusted if they were
			return False
	else:
		#is a duplicate
		return True

def tidy_PAIRED_dict(dict,current_position1,current_position2,read_length):
	#We'll use the following logic to remove
	#reads that fall out of the range that could possibly be
	#duplicates we'll remove purge from the dictionary
	new = {k:v for k,v in dict.items() if k[4] >= (int(current_position2) + read_length) or k[3] >= (int(current_position1) - read_length) }		
	return new
	
######################
# GLOBAL VARIABLES 2 #
######################

max_read_length = 0	# maximum read length
withoutPOS = 0			# unmapped read counter
total_reads = 0			# read counter
ForwardStrands = dict() # dictionaries for forward strand
ReverseStrands = dict() # dictionaries for reverse strand

#############################
#  		MAIN DEDUPER        #
#############################
def Singlend():
	global max_read_length
	global total_reads
	global withoutPOS
	global ForwardStrands
	global ReverseStrands
	pure_read_count = 0
	with open(file, "r") as fh, open(file+"_deduped","w") as wh:
		
		for i,line in enumerate(fh):	
			if "@" != line.split()[0][0]:
				total_reads += 1
				lineT = line.split("\t") #split by tabs
				#[0 1 2 3 ...       ] strings
				flag = int(lineT[1]) #binary flag, stranded/paired or single end
				UMI = lineT[0].split(":")[-1] #unique molecular identifier
				CIGAR = lineT[5] #soft clipping, other stuff
				POS = lineT[3]   #position
				max_read_length = max(len(lineT[9]),max_read_length)
				tidy_it = (total_reads % max_read_length == 0) #tidy every so often to keep dictionaries tidy
				CHRM = lineT[2]
				if POS.isdigit() == False:
					#we got trouble here, not mapped, can't dedupe
					#abort or skip?
					withoutPOS	+= 1
					continue #skip, which is equivalent to being a duplicate
				
				if (UMI_list_flag): #UMI list provided
					if (UMI not in umi_list):
						#try to repair UMI?
						if strict > 0:					
							try:
								UMI_correction = umi_list[sorted([(hamming_distance(x,UMI),i) for i,x in enumerate(umi_list) if hamming_distance(x,UMI) <= strict])[0][1]]
								line = line.rstrip() + "\tCU:Z:"+UMI_correction+"\n" # add a tag with the corrected UMI
								UMI = UMI_correction
							except:
								continue #couldn't find correction
						else:
							#not in list and no correction (skip)
							continue
						
					#UMI is in the list!
				#else no UMI list
				
				
				
				if ((flag & 16) == 16): #check flag for reverse compliment
					#sequence reverse complimented 
					#hard mode
					X = CIGAR_REV_PARSE(CIGAR)
					Y = X + int(POS) #this is our reference which will go along with our dictionary to keep track when we can get rid of the entries and so on
					
					if not is_Duplicate(ReverseStrands,UMI,CHRM,Y):
						#it is NOT a duplicate
						add_read(ReverseStrands,line,Y,UMI,CHRM)
						wh.write(line)
						pure_read_count += 1
						#finally we'll check every once in a 
						#while if we can remove some
						#of the reads in the dictionary
						if tidy_it:
							ReverseStrands = tidy_REV_dict(ReverseStrands,POS,max_read_length)
						#else it is a duplicate and we do nothing
				
				else:
					#forward strand
					#easy mode
					X = CIGAR_FOR_PARSE(CIGAR) #get soft clipping adjustment
					Y = int(POS) - X #this is our reference which will go along with our dictionary to keep track when we can get rid of the entries and so on	
					if not is_Duplicate(ForwardStrands,UMI,CHRM,Y):
						#is not a duplicate
						#add to ForwardStrands and write to file
						add_read(ForwardStrands,line,Y,UMI,CHRM)
						wh.write(line)
						pure_read_count += 1
						#finally we'll check every once in a 
						#while if we can remove some
						#of the reads in the dictionary
						if tidy_it:
							ForwardStrands = tidy_FOR_dict(ForwardStrands,POS,max_read_length)
						#else it is a duplicate and we do nothing
			else:
				wh.write(line) #write header lines
	#print(pure_read_count)
	return None#done

######################
# GLOBAL VARIABLES 3 #
######################

paired_Dict = dict() 



	
def Paired():

	"""
	check 9:78,175,324-78,175,364 in test file
	"""
	global max_read_length
	global total_reads
	global paired_Dict
	skip_lines = []
	total_lines = 0
	forwardAdjust = 0
	reverseAdjust = 0
	firstIn = False
	first_UMI_fail = False
	max_search = 0
	total = 0
	FLAG = 0
	dup_pair = []
	pair_not_found = []
	cur_pos_max = 0
	with open(file, "r") as fh, open(file+"_deduped","w") as wh:
		previous = ""
		cur_count = 0
		#get some information about searching for the 
		#next pair
		for i,line in enumerate(fh):
			flag = line.split("\t")[1]
			if (flag != previous):
				cur_count = 1
				previous = flag
				
			else:
				cur_count += 1
				max_search = max(cur_count,max_search)	
			total = i+1
		max_search *= 2
		fh.seek(0)
		while True:
			line1 = fh.readline()
			
			total_lines += 1
			if total_lines > total:
				break
			if total_lines in skip_lines:
				skip_lines.pop(0)
				continue
			if "@" in line1.split()[0][0]:
				wh.write(line1)
			else:
				#get info
				total_reads += 1
				line1T = line1.split("\t") #split by tabs
				#[0 1 2 3 ...       ] strings
				flag1 = int(line1T[1]) #binary flag, stranded/paired or single end
				UMI = line1T[0].split(":")[-1] #unique molecular identifier
				UMI1 = UMI.split("^")[0]
				UMI2 = UMI.split("^")[1]
				CIGAR1 = line1T[5] #soft clipping, other stuff
				PNEXT = line1T[7]
				POS1 = line1T[3]   #position
				max_read_length = max(len(line1T[9]),max_read_length)
				CHRM1 = line1T[2]
				tidy_it = (total_reads % 100) == 0
				PairedRead = ""
				first_UMI_fail = False
				firstIn = ((flag1 & 64) == 64)
				if (UMI_list_flag): #UMI list provided
					UMI1c,pass1,check1 = UmiListCorrection(UMI1)
					UMI2c,pass2,check2 = UmiListCorrection(UMI2)
					if not (pass1 and pass2):
						#must find pair and kill it
						first_UMI_fail = True
					elif check1 or check2:
						#umi corrected, add a TAG
						line1 = line1.rstrip() + "\tCU:Z:"+UMI1c+"^"+UMI2c+"\n" # add a tag with the corrected UMI
					UMI1 = UMI1c
					UMI2 = UMI2c	
				#find paired read
				remember = fh.tell() #remember where we were
				
				if (firstIn):
					for i in range(max_search):
						line = fh.readline()
						if line == "":
							break
						#get info
						line2T = line.split("\t") #split by tabs
						#[0 1 2 3 ...       ] strings
						flag2 = int(line2T[1]) #binary flag, stranded/paired or single end
						#find 2nd in
						if ((flag2 & 128) == 128) and (PNEXT == line2T[3]):
							UMIs = line2T[0].split(":")[-1] #unique molecular identifier
							UMI1s = UMIs.split("^")[0]
							UMI2s = UMIs.split("^")[1]
							if UMI1s == UMI1 and UMI2s == UMI2:
								#pair found
								
								PairedRead = line
								skip_lines.append(total_lines+i+1)
								break #break out of finding pair
							elif UMI_list_flag:
								pass1 = (hamming_distance(UMI1,UMI1s) <= strict)
								pass2 = (hamming_distance(UMI2,UMI2s) <= strict)
								if not (pass1 and pass2):
									continue #failed UMI test, find another
								else:
									PairedRead = line
									skip_lines.append(total_lines+i+1)
									break
							elif first_UMI_fail:
								if (hamming_distance(UMI1,UMI1s) <= 2 and hamming_distance(UMI2,UMI1s) <= 2):
									skip_lines.append(total_lines+i+1)
									break
							elif i == (max_search-1):
								#couldn't find a pair
								pair_not_found.append((CHRM1,POS1))
					
					
				else:
					for i in range(max_search):
						line = fh.readline()
						if line == "":
							break
						#get info
						line2T = line.split("\t") #split by tabs
						#[0 1 2 3 ...       ] strings
						flag2 = int(line2T[1]) #binary flag, stranded/paired or single end
						
						#find 2nd in
						if ((flag2 & 64) == 64) and (PNEXT == line2T[3]):
							UMIs = line2T[0].split(":")[-1] #unique molecular identifier
							UMI1s = UMIs.split("^")[0]
							UMI2s = UMIs.split("^")[1]
							if UMI1s == UMI1 and UMI2s == UMI2:
								#pair found
								
								PairedRead = line
								skip_lines.append(total_lines+i+1)
								break #break out of finding pair
							elif UMI_list_flag:
								pass1 = (hamming_distance(UMI1,UMI1s) <= strict)
								pass2 = (hamming_distance(UMI2,UMI2s) <= strict)
								if not (pass1 and pass2):
									continue #failed UMI test, find another
								else:
									PairedRead = line
									skip_lines.append(total_lines+i+1)
									break
							elif first_UMI_fail:
								if (hamming_distance(UMI1,UMI1s) <= 2 and hamming_distance(UMI2,UMI1s) <= 2):
									skip_lines.append(total_lines+i+1)
									break
							elif i == (max_search-1):
								pair_not_found.append((CHRM1,POS1))
						
					
				fh.seek(remember)
				if PairedRead == "" or first_UMI_fail:
					continue # couldn't find second read or first UMIs not in list
				#get info of 2nd read
				line2T = PairedRead.split("\t") #split by tabs
				#[0 1 2 3 ...       ] strings
				flag2 = int(line2T[1]) #binary flag, stranded/paired or single end
				CIGAR2 = line2T[5] #soft clipping, other stuff
				POS2 = line2T[3]   #position
				max_read_length = max(len(line1T[9]),max_read_length,len(line2T[9]))
				CHRM2 = line2T[2]
				firstReverse = ((flag1 & 16) == 16) #sequence1 is reverse complimented
				
				if firstReverse:
					reverseAdjust = int(POS1) + CIGAR_REV_PARSE(CIGAR1)
					forwardAdjust = int(POS2) - CIGAR_FOR_PARSE(CIGAR2)
				else:
					reverseAdjust = int(POS2) + CIGAR_REV_PARSE(CIGAR2)
					forwardAdjust = int(POS1) - CIGAR_FOR_PARSE(CIGAR1)
				if firstIn:
					FLAG = flag1
				else:
					FLAG = flag2
				if not isPairedDuplicate(paired_Dict,UMI1,UMI2,CHRM1,forwardAdjust,reverseAdjust,FLAG):
					add_paired_read(paired_Dict,PairedRead,line1,forwardAdjust,reverseAdjust,UMI1,UMI2,CHRM1,FLAG)
					wh.write(line1) #keeps first
					wh.write(PairedRead)
				
				
				if tidy_it:
					paired_Dict = tidy_PAIRED_dict(paired_Dict,POS1,POS2,max_read_length)
					
					
					#print("Percent Complete:", float(total_lines)/float(total), total_lines,"/",total)
				
	
	return None#done
	
	
if not paired:
	Singlend()
else:
	Paired()