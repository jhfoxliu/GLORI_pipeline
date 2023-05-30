#!bin/usr/env python

#Jianheng Liu @ Zhanglab, SYSU
#Feb, 2018
#Email: liujh26@mail2.sysu.edu.cn
#Usage: Map reads to HISAT2 indexed genome
#Input: [.fastq] [HISAT2 index]

__doc__ = """
Usage: Map reads to HISAT2 indexed genome.
Note that the input should be stranded. If the conversion
rates are very low, this program might take too much RAM.
In such a situation, please split the fastq input and merge
BAMs later.
"""


import os,sys,time,signal,argparse,pysam,multiprocessing
from Bio import SeqIO
from Bio.Seq import reverse_complement
from time import gmtime, strftime
from collections import defaultdict
from pysam import qualities_to_qualitystring
import subprocess


def read_parameters(para_file):
	parameters = {}
	if para_file:
		with open(para_file, 'r') as input:
			for line in input.readlines():
				if line.startswith("#") == False:
					line = line.split(" ")
					parameters[line[0]] = line[1]
	return parameters

def fastq_converter_1(fin,fout,method):
	sys.stderr.write("[%s]Converting %s, %s...\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()),fin,method))
	total_reads = 0
	output_dict = {}
	if method == "A2G":
		original = "A"
		convert = "G"
	elif method == "T2C":
		original = "T"
		convert = "C"
	with open(fin,'r') as input, open(fout,'w') as output:
		line = input.readline()
		i = 0
		while (line):
			if i == 0:
				output.write(line)
				name = line.strip("@").strip("\n").split( )[0]
				i += 1
			elif i == 1:
				line = line.strip()
				if original in line:
					output_dict[name] = line
					output.write(line.replace(original,convert))
					output.write("\n")
				else:
					output.write(line)
					output.write("\n")
				i += 1
			elif i == 2:
				output.write(line)
				i +=1
			elif i == 3:
				output.write(line)
				i = 0
				total_reads += 1
			line = input.readline()
	return output_dict,total_reads
	
def read_m6A_position(fwd,rev=None):
	total_reads = 0
	fwd_read_dict = {}
	rev_read_dict = {}
	if rev is not None:
		with open(fwd,'r') as fastq_fwd, open(rev,'r') as fastq_rev:
			line_fwd = fastq_fwd.readline()
			line_rev = fastq_rev.readline()
			n = 0
			while (line_fwd):
				if n == 0:
					seq_name = line_fwd.strip("@").strip("\n").split( )[0]
					n=n+1
					total_reads += 1
				elif n == 1:
					fwd_C = line_fwd.strip("\n").count("A")
					rev_G = line_rev.strip("\n").count("T")
					if fwd_C != 0:
						fwd_read_dict[seq_name] = line_fwd.strip("\n")
					if rev_G != 0:
						rev_read_dict[seq_name] = line_rev.strip("\n")
					n=n+1
				elif n == 2:
					n=n+1
				elif n == 3:
					n=0
				line_fwd = fastq_fwd.readline()
				line_rev = fastq_rev.readline()
	else:
		with open(fwd,'r') as fastq_fwd:
			line_fwd = fastq_fwd.readline()
			n = 0
			while (line_fwd):
				if n == 0:
					seq_name = line_fwd.strip("@").strip("\n").split( )[0]
					n=n+1
					total_reads += 1
				elif n == 1:
					fwd_C = line_fwd.strip("\n").count("A")
					if fwd_C != 0:
						fwd_read_dict[seq_name] = line_fwd.strip("\n")
					n=n+1
				elif n == 2:
					n=n+1
				elif n == 3:
					n=0
				line_fwd = fastq_fwd.readline()
	return fwd_read_dict,rev_read_dict,total_reads
	
def mapping_string(parameters,fwd,rev):
	if "-p" not in parameters:
		parameters["-p"] = '5'
	if "-k" not in parameters:
		parameters["-k"] = '10'
	if "--reorder" not in parameters:
		parameters["--reorder"] = None
	if fwd and rev:
		parameters["-1"] = fwd
		parameters["-2"] = rev
		if "--rna-strandness" not in parameters:
			parameters["--rna-strandness"] = "FR"
		if "--no-mixed" not in parameters:
			parameters["--no-mixed"] = None
		if "--fr" not in parameters:
			parameters["--fr"] = None
	elif fwd and not rev:
		parameters["-U"] = fwd
		if "--rna-strandness" not in parameters:
			parameters["--rna-strandness"] = "F"
	if "--no-temp-splicesite" not in parameters:
		parameters["--no-temp-splicesite"] = None
		
	# if options.avoid_psedu and "--avoid-pseudogene" not in parameters:
		# parameters["--avoid-pseudogene"] = None
		
	para = []
	for key,values in parameters.items():
		if values is not None:
			para.extend([key,values])
		else:
			para.extend([key])
	
	para = "|".join(para)
	return para
	
def mapping(index,method,hisat2_path,parameters,preifx):
	if parameters:
		parameters = parameters.split("|")
	else:
		parameters = []
	summary_name = method + "hisat2.summary"

	if method == "A2G":
		sam = preifx + ".A2G.sam"
	elif method == "T2C":
		sam = preifx + ".T2C.sam"
	cmd = [hisat2_path+"/hisat2", \
		    "-x",index,\
			"-S",sam] + \
		   parameters

	hisat2_out = subprocess.Popen(cmd,\
								  stdin=subprocess.PIPE,\
								  stdout=subprocess.PIPE,\
								  stderr=subprocess.PIPE,\
								  shell=False)
	stdout,stderr = hisat2_out.communicate()
	# print stdout
	sys.stderr.write("[%s]%s report:\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()),method))
	sys.stderr.write(stderr)
	sys.stderr.write("\n")
	#stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=False
	# hisat2.wait()
def set_tag(segment_tag,genome):
	segment_tag.set_tag("YG",genome)
	if segment_tag.is_read1:
		segment_tag.set_tag("YR","A2G")
	elif segment_tag.is_read2:
		segment_tag.set_tag("YR","T2C")
	else:
		if segment_tag.is_reverse:
			segment_tag.set_tag("YR","T2C")
		else:
			segment_tag.set_tag("YR","A2G")
	return segment_tag
	
def check_mapping(segment_check):
	if segment_check.get_tag("YG") == "A2G" and segment_check.get_tag("XS") == "+" and segment_check.get_tag("AS") >= options.min_AS:
		if segment_check.is_read1 and segment_check.is_reverse == False:
			return True,segment_check
		elif segment_check.is_read2 and segment_check.is_reverse == True:
			return True,segment_check
		elif segment_check.is_paired == False and segment_check.is_reverse == False:
			return True,segment_check
		else:
			return False,segment_check
	elif segment_check.get_tag("YG") == "T2C" and segment_check.get_tag("XS") == "-" and segment_check.get_tag("AS") >= options.min_AS:
		if segment_check.is_read1 and segment_check.is_reverse == True:
			return True,segment_check
		elif segment_check.is_read2 and segment_check.is_reverse == False:
			return True,segment_check
		elif segment_check.is_paired == False and segment_check.is_reverse == True:
			return True,segment_check
		else:
			return False,segment_check
	else:
		return False,segment_check
		
def check_list(TMP,genome,discordant_A2G,discordant_T2C):
	passed = []
	for segment in TMP:
		segment = set_tag(segment,genome)
		status, segment = check_mapping(segment)
		if status == True:
			passed.append(segment)
		else:
			if segment.get_tag("YG") == "A2G":
				discordant_A2G.append(segment)
			elif segment.get_tag("YG") == "T2C":
				discordant_T2C.append(segment)
	return passed,discordant_A2G,discordant_T2C
	
def to_unmapped(segment,unal_read1,unal_read2):
	if segment.is_read1:
		if segment.query_name in fwd_read_dict:
			seq = fwd_read_dict[segment.query_name]
			qual = qualities_to_qualitystring(segment.query_qualities)
			if segment.is_reverse:
				qual = qual[::-1]
		else:
			seq = segment.query_sequence
			qual = qualities_to_qualitystring(segment.query_qualities)
			if segment.is_reverse:
				seq = reverse_complement(seq)
				qual = qual[::-1]
		unal_read1.write("".join(["@",segment.query_name,"\n",seq,'\n+\n',qual,'\n']))
	elif segment.is_read2:
		if segment.query_name in rev_read_dict:
			seq = rev_read_dict[segment.query_name]
			qual = qualities_to_qualitystring(segment.query_qualities)
			if segment.is_reverse:
				qual = qual[::-1]
		else:
			seq = segment.query_sequence
			qual = qualities_to_qualitystring(segment.query_qualities)
			if segment.is_reverse:
				seq = reverse_complement(seq)
				qual = qual[::-1]
		unal_read2.write("".join(["@",segment.query_name,"\n",seq,'\n+\n',qual,'\n']))
	else: #single end and fully unmapped
		if segment.query_name in fwd_read_dict:
			seq = fwd_read_dict[segment.query_name]
			qual = qualities_to_qualitystring(segment.query_qualities)
			if segment.is_reverse:
				qual = qual[::-1]
		else:
			seq = segment.query_sequence
			qual = qualities_to_qualitystring(segment.query_qualities)
			if segment.is_reverse:
				seq = reverse_complement(seq)
				qual = qual[::-1]
		unal_read1.write("".join(["@",segment.query_name,"\n",seq,'\n+\n',qual,'\n']))

def to_unmapped_discordant(discordant,unal_read1,unal_read2):
	read1 = None
	read2 = None
	discordant_iter = iter(discordant)
	while read1 is None or read2 is None:
		segment = next(discordant_iter)
		if not read1 and not read2:
			if segment.is_read1:
				to_unmapped(segment,unal_read1,unal_read2)
				read1 = 1
			elif segment.is_read2:
				to_unmapped(segment,unal_read1,unal_read2)
				read2 = 1
			else:
				to_unmapped(segment,unal_read1,unal_read2)
				read1 = 1
				read2 = 1
		elif not read1 and read2:
			if segment.is_read1:
				to_unmapped(segment,unal_read1,unal_read2)
				read1 = 1
		elif not read2 and read1:
			if segment.is_read2:
				to_unmapped(segment,unal_read1,unal_read2)
				read2 = 1
		else:
			break
	
def to_bam_genome(segments,output):
	for segment in segments:
		if segment.is_read1:
			seq = fwd_read_dict.get(segment.query_name)
		elif segment.is_read2:
			seq = rev_read_dict.get(segment.query_name)
		else:
			seq = fwd_read_dict.get(segment.query_name)
		if seq is not None:
			if segment.is_reverse:
				seq = reverse_complement(str(seq))
			qual = segment.query_qualities
			mpq = segment.mapping_quality
			segment_output = segment
			segment_output.query_sequence = seq
			segment_output.mapping_quality = mpq
			segment_output.query_qualities = qual
			output.write(segment_output)
		else:
			output.write(segment)
	
def check_pairs(has_read1,has_read2,segment):
	if segment.is_read1:
		has_read1 = True
	elif segment.is_read2:
		has_read2 = True
	return has_read1,has_read2
	
def validate_read_pairs(a,b):
	if (a.is_read1 and b.is_read2) or (a.is_read2 and b.is_read1):
		if a.reference_name == b.reference_name:
			return True
		else:
			return False
	else:
		return False

def remove_high_poly(segment):
	A = 0.0
	T = 0.0
	C = 0.0
	G = 0.0
	for i in segment.query_sequence:
		if i == "A":
			A += 1
		elif i == "T":
			T += 1
		elif i == "C":
			C += 1
		elif i == "G":
			G += 1
	Total = A + T + C +G
	if A/Total >= 0.8 or T/Total >= 0.8 or C/Total >= 0.8 or G/Total >= 0.8:
		return True
	else:
		return False


def mapping_result(unique_A2G_tmp,unique_T2C_tmp,unmapped_A2G,unmapped_T2C,multi_A2G_tmp,multi_T2C_tmp,\
                   output,unal_read1,unal_read2,multi_output,mode):
	global multimapper_num,unmapped_num,A2G_num,T2C_num
	discordant_A2G,discordant_T2C = [],[]
	unique_A2G,discordant_A2G,discordant_T2C = check_list(unique_A2G_tmp,"A2G",discordant_A2G,discordant_T2C)
	unique_T2C,discordant_A2G,discordant_T2C = check_list(unique_T2C_tmp,"T2C",discordant_A2G,discordant_T2C)
	multi_A2G,discordant_A2G,discordant_T2C = check_list(multi_A2G_tmp,"A2G",discordant_A2G,discordant_T2C)
	multi_T2C,discordant_A2G,discordant_T2C = check_list(multi_T2C_tmp,"T2C",discordant_A2G,discordant_T2C)
	
	if mode == "PE": #no more multimapping
		if len(unique_A2G) == 1:
			discordant_A2G.extend(unique_A2G)
			unique_A2G = []
		if len(unique_T2C) == 1:
			discordant_T2C.extend(unique_T2C)
			unique_T2C = []
		if len(multi_A2G) == 2:
			has_read1,has_read2 = False,False
			# has_read1,has_read2 = check_pairs(has_read1,has_read2,multi_A2G[0])
			#has_read1,has_read2 = check_pairs(has_read1,has_read2,multi_A2G[1])
			if validate_read_pairs(multi_A2G[0],multi_A2G[1]) == True:#has_read1 and has_read2:
				while multi_A2G:
					segment = multi_A2G.pop()
					segment.is_secondary = False
					segment.set_tag("NH",1)
					unique_A2G.append(segment)
			else:
				discordant_A2G.extend(multi_A2G)
				multi_A2G = []
		if len(multi_T2C) == 2:
			#has_read1,has_read2 = False,False
			#has_read1,has_read2 = check_pairs(has_read1,has_read2,multi_T2C[0])
			#has_read1,has_read2 = check_pairs(has_read1,has_read2,multi_T2C[1])
			if validate_read_pairs(multi_T2C[0],multi_T2C[1]) == True:#has_read1 and has_read2:
				while multi_T2C:
					segment = multi_T2C.pop()
					segment.is_secondary = False
					segment.set_tag("NH",1)
					unique_T2C.append(segment)
			else:
				discordant_T2C.extend(multi_T2C)
				multi_T2C = []
	elif mode == "SE":
		if len(multi_A2G) == 1:
			while multi_A2G:
				segment = multi_A2G.pop()
				segment.is_secondary = False
				segment.set_tag("NH",1)
				unique_A2G.append(segment)
		if len(multi_T2C) == 1:
			while multi_T2C:
				segment = multi_T2C.pop()
				segment.is_secondary = False
				segment.set_tag("NH",1)
				unique_T2C.append(segment)

	if unmapped_A2G and unmapped_T2C:
		for segment in unmapped_A2G:
			to_unmapped(segment,unal_read1,unal_read2)
		unmapped_num += 1
	elif multi_A2G or multi_T2C:
		to_bam_genome(multi_A2G,multi_output)
		to_bam_genome(multi_T2C,multi_output)
		to_bam_genome(unique_A2G,multi_output)
		to_bam_genome(unique_T2C,multi_output)
		if options.no_multi_fastq == False:
			if multi_A2G:
				to_unmapped_discordant(multi_A2G+discordant_A2G,unal_read1,unal_read2)
			elif multi_T2C:
				to_unmapped_discordant(multi_T2C+discordant_T2C,unal_read1,unal_read2)
		multimapper_num += 1
	elif unique_A2G and unique_T2C:
		''' experimental function '''
		if options.fully_unique == False:
			if mode == "PE":
				if unique_A2G[0].is_proper_pair == True and unique_T2C[0].is_proper_pair == False:
					for i in unique_A2G:
						is_poly = remove_high_poly(i)
						if is_poly == True:
							to_unmapped_discordant(unique_A2G,unal_read1,unal_read2)
							unmapped_num += 1
							return None
					to_bam_genome(unique_A2G,output)
					A2G_num += 1
				elif unique_A2G[0].is_proper_pair == False and unique_T2C[0].is_proper_pair == True:
					for i in unique_T2C:
						is_poly = remove_high_poly(i)
						if is_poly == True:
							to_unmapped_discordant(unique_T2C,unal_read1,unal_read2)
							unmapped_num += 1
							return None
					to_bam_genome(unique_T2C,output)
					T2C_num += 1
				else:
					A2G_AS = unique_A2G[0].get_tag("AS") + unique_A2G[1].get_tag("AS")
					T2C_AS = unique_T2C[0].get_tag("AS") + unique_T2C[1].get_tag("AS")
					if A2G_AS - T2C_AS > options.diff_AS:
						for i in unique_A2G:
							is_poly = remove_high_poly(i)
							if is_poly == True:
								to_unmapped_discordant(unique_A2G,unal_read1,unal_read2)
								unmapped_num += 1
								return None
						to_bam_genome(unique_A2G,output)
						A2G_num +=1 
					elif T2C_AS - A2G_AS > options.diff_AS:
						for i in unique_T2C:
							is_poly = remove_high_poly(i)
							if is_poly == True:
								to_unmapped_discordant(unique_T2C,unal_read1,unal_read2)
								unmapped_num += 1
								return None
						to_bam_genome(unique_T2C,output)
						T2C_num +=1 
					else:
						to_bam_genome(unique_A2G,multi_output)
						to_bam_genome(unique_T2C,multi_output)
						if options.no_multi_fastq == False:
							to_unmapped_discordant(unique_A2G,unal_read1,unal_read2)
						multimapper_num += 1
			elif mode == "SE":
				A2G_AS = unique_A2G[0].get_tag("AS")
				T2C_AS = unique_T2C[0].get_tag("AS")
				if A2G_AS - T2C_AS > options.diff_AS:
					for i in unique_A2G:
						is_poly = remove_high_poly(i)
						if is_poly == True:
							to_unmapped_discordant(unique_A2G,unal_read1,unal_read2)
							unmapped_num += 1
							return None
					to_bam_genome(unique_A2G,output)
					A2G_num +=1 
				elif A2G_AS - T2C_AS > options.diff_AS:
					for i in unique_T2C:
						is_poly = remove_high_poly(i)
						if is_poly == True:
							to_unmapped_discordant(unique_T2C,unal_read1,unal_read2)
							unmapped_num += 1
							return None
					to_bam_genome(unique_T2C,output)
					T2C_num +=1 
				else:
					to_bam_genome(unique_A2G,multi_output)
					to_bam_genome(unique_T2C,multi_output)
					if options.no_multi_fastq == False:
						to_unmapped_discordant(unique_A2G,unal_read1,unal_read2)
					multimapper_num += 1
		else:
			to_bam_genome(unique_A2G,multi_output)
			to_bam_genome(unique_T2C,multi_output)
			if options.no_multi_fastq == False:
				to_unmapped_discordant(unique_A2G,unal_read1,unal_read2)
			multimapper_num += 1
	elif unique_A2G and not unique_T2C:
		for i in unique_A2G:
			is_poly = remove_high_poly(i)
			if is_poly == True:
				to_unmapped_discordant(unique_A2G,unal_read1,unal_read2)
				unmapped_num += 1
				return None
		to_bam_genome(unique_A2G,output)
		A2G_num += 1
	elif unique_T2C and not unique_A2G:
		for i in unique_T2C:
			is_poly = remove_high_poly(i)
			if is_poly == True:
				to_unmapped_discordant(unique_T2C,unal_read1,unal_read2)
				unmapped_num += 1
				return None
		to_bam_genome(unique_T2C,output)
		T2C_num += 1
	elif any([unique_A2G,unique_T2C,multi_A2G,multi_T2C]) == False:
		if discordant_A2G:
			to_unmapped_discordant(discordant_A2G,unal_read1,unal_read2)
		elif discordant_T2C:
			to_unmapped_discordant(discordant_T2C,unal_read1,unal_read2)
		unmapped_num += 1

def itering_groups(A2G_iter,T2C_iter,A2G_next,T2C_next,\
				   output,unal_read1,unal_read2,multi_output,mode):
		unique_A2G = []
		unique_T2C = []
		unmapped_A2G = []
		unmapped_T2C = []
		multi_A2G = []
		multi_T2C = []
		#first read
		try:
			if not A2G_next and not T2C_next:
				A2G_read = A2G_iter.__next__()
				T2C_read = T2C_iter.__next__()
			else:
				A2G_read = A2G_next
				T2C_read = T2C_next
		except StopIteration:
			return None,None,None,None
			
		last_read = A2G_read.query_name
		A2G_group = A2G_read.get_tag("NH") if A2G_read.has_tag("NH") else None

		try:
			while A2G_read.query_name == last_read:
				if A2G_group == 1: #A2G unique
					unique_A2G.append(A2G_read)
				elif A2G_group is None: #Unmapped
					unmapped_A2G.append(A2G_read)
				else: #multimapper:
					multi_A2G.append(A2G_read)
				A2G_read = A2G_iter.__next__()
			A2G_next = A2G_read
		except StopIteration:
			pass
			
		T2C_group = T2C_read.get_tag("NH") if T2C_read.has_tag("NH") else None
		try:
			while T2C_read.query_name == last_read:
				if T2C_group == 1: #T2C unique
					unique_T2C.append(T2C_read)
				elif T2C_group is None: #Unmapped
					unmapped_T2C.append(T2C_read)
				else: #multimapper:
					multi_T2C.append(T2C_read)
				T2C_read = T2C_iter.__next__()
			T2C_next = T2C_read
			mapping_result(unique_A2G,unique_T2C,unmapped_A2G,unmapped_T2C,multi_A2G,multi_T2C,\
						   output,unal_read1,unal_read2,multi_output,mode)
			return A2G_iter,T2C_iter,A2G_next,T2C_next
		except StopIteration:
			mapping_result(unique_A2G,unique_T2C,unmapped_A2G,unmapped_T2C,multi_A2G,multi_T2C,\
						   output,unal_read1,unal_read2,multi_output,mode)
			A2G_next = None
			T2C_next = None
			return A2G_iter,T2C_iter,A2G_next,T2C_next

def sam_handler_PE():
	unal_1 = options.fwd.replace(".fastq",".unmapped.fastq")
	unal_2 = options.rev.replace(".fastq",".unmapped.fastq")
	if not options.continue_prefix:
		hisat2_prefix =  "hisat2_" + str(os.getpid())
		hisat2_A2G_SAM = hisat2_prefix + ".A2G.sam"
		hisat2_T2C_SAM = hisat2_prefix + ".T2C.sam"
	else:
		hisat2_A2G_SAM = options.continue_prefix+".A2G.sam"
		hisat2_T2C_SAM = options.continue_prefix+".T2C.sam"
		
	# with pysam.AlignmentFile("hisat2.A2G.sam","r") as A2G_SAM, \
		 # pysam.AlignmentFile("hisat2.T2C.sam","r") as T2C_SAM, \
		 # pysam.AlignmentFile(options.output+".bam","wb",template=A2G_SAM) as output, \
		 # open(unal_1,'w') as unal_read1, \
		 # open(unal_2,'w') as unal_read2, \
		 # pysam.AlignmentFile(options.output+".multimappers.bam","wb",template=A2G_SAM) as multi_output:
	with pysam.AlignmentFile(hisat2_A2G_SAM,"r") as A2G_SAM, \
		 pysam.AlignmentFile(hisat2_T2C_SAM,"r") as T2C_SAM, \
		 pysam.AlignmentFile(options.output+".bam","wb",template=A2G_SAM) as output, \
		 open(unal_1,'w') as unal_read1, \
		 open(unal_2,'w') as unal_read2, \
		 pysam.AlignmentFile(options.output+".multimappers.bam","wb",template=A2G_SAM) as multi_output:
		A2G_iter = A2G_SAM.fetch(until_eof=True)
		T2C_iter = T2C_SAM.fetch(until_eof=True)
		A2G_next = None
		T2C_next = None
		while (A2G_iter and T2C_iter):
			A2G_iter,T2C_iter,A2G_next,T2C_next = itering_groups(A2G_iter,T2C_iter,A2G_next,T2C_next,output,unal_read1,unal_read2,multi_output,"PE")

def itering_groups_SE(A2G_iter,A2G_next,\
				      output,unal_read1,multi_output,mode):
		unique_A2G = []
		unique_T2C = []
		unmapped_A2G = []
		unmapped_T2C = []
		multi_A2G = []
		multi_T2C = []
		#first read
		try:
			if not A2G_next and not T2C_next:
				A2G_read = A2G_iter.__next__()
			else:
				A2G_read = A2G_next
		except StopIteration:
			return None,None

		last_read = A2G_read.query_name
		A2G_group = A2G_read.get_tag("NH") if A2G_read.has_tag("NH") else None
		
		try:
			while A2G_read.query_name == last_read:
				if A2G_group == 1: #A2G unique
					unique_A2G.append(A2G_read)
				elif A2G_group is None: #Unmapped
					unmapped_A2G.append(A2G_read)
				else: #multimapper:
					multi_A2G.append(A2G_read)
				A2G_read = A2G_iter.__next__()
			A2G_next = A2G_read
		except StopIteration:
			mapping_result(unique_A2G,unique_T2C,unmapped_A2G,unmapped_T2C,multi_A2G,multi_T2C,\
						   output,unal_read1,unal_read2,multi_output,mode)
			A2G_next = None
			return A2G_iter,A2G_next
			
def sam_handler_SE():
	unal_1 = options.fwd.replace(".fastq",".unmapped.fastq")
	unal_2 = "single_end.ignore"
	if not options.continue_prefix:
		hisat2_prefix = "hisat2_" + str(os.getpid())
		hisat2_A2G_SAM = hisat2_prefix + ".A2G.sam"
		hisat2_T2C_SAM = hisat2_prefix + ".T2C.sam"
	else:
		hisat2_A2G_SAM = options.continue_prefix+".A2G.sam"
		hisat2_T2C_SAM = options.continue_prefix+".T2C.sam"
	with pysam.AlignmentFile(hisat2_A2G_SAM,"r") as A2G_SAM, \
		 pysam.AlignmentFile(hisat2_T2C_SAM,"r") as T2C_SAM, \
		 pysam.AlignmentFile(options.output+".bam","wb",template=A2G_SAM) as output, \
		 open(unal_1,'w') as unal_read1, \
		 open(unal_2,'w') as unal_read2, \
		 pysam.AlignmentFile(options.output+".multimappers.bam","wb",template=A2G_SAM) as multi_output:
		A2G_iter = A2G_SAM.fetch(until_eof=True)
		T2C_iter = T2C_SAM.fetch(until_eof=True)
		A2G_next = None
		T2C_next = None
		while (A2G_iter and T2C_iter):
			A2G_iter,T2C_iter,A2G_next,T2C_next = itering_groups(A2G_iter,T2C_iter,A2G_next,T2C_next,output,unal_read1,unal_read2,multi_output,"SE")
	os.remove(unal_2)
	
def signal_handler(sig,frame):
	try:
		pool.terminate()
	except NameError:
		pass
	sys.exit()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog="BS_hisat2",fromfile_prefix_chars='@',description=__doc__,formatter_class=argparse.RawTextHelpFormatter)
	#Required
	group_required = parser.add_argument_group("Required")
	group_required.add_argument("--fwd","-F",dest="fwd",required=True,help="Forward read (A less read)")
	group_required.add_argument("--rev","-R",dest="rev",required=False,help="Reverse read (T less read)")
	group_required.add_argument("--index","-I",dest="index",required=True,help="Path to index folder")
	group_required.add_argument("--output","-o",dest="output",required=True,help="Output prefix")
	#Hisat2
	group_hisat2 = parser.add_argument_group("Hisat2")
	group_hisat2.add_argument("--index-prefix",dest="index_prefix",required=False,default="HISAT2",help="Index will looks like /path/A2G/[prefix]_A2G and /path/T2C/[prefix]_T2C")
	group_hisat2.add_argument("--hisat2-path",dest="hisat2_path",help="Path to hisat2 bin")
	group_hisat2.add_argument("--hisat2-param",dest="hisat2_param",help="Index building parameters")
	#Filter
	group_filter = parser.add_argument_group("Filter")
	group_filter.add_argument("--min-AS",dest="min_AS",type=int,default=-10,help="minimal AS")
	group_filter.add_argument("--diff-AS",dest="diff_AS",type=int,default=9,help="AS different between best and secondary mapping")
	group_filter.add_argument("--fully-unique",dest="fully_unique",default=False,action="store_true",help="Reads should be aligned to A2G or T2C by only ONE time")
	#Steps
	group_steps = parser.add_argument_group("Steps")
	group_steps.add_argument("--continue-prefix",dest="continue_prefix",help="If given, process the sam files: [preifx].A2G.sam and [prefix].T2C.sam. If not given, prefix=hisat2_[pid]")
	group_steps.add_argument("--skip-convert",dest="skip_convert",default=False,action="store_true",help="Skip converting fastq")
	group_steps.add_argument("--skip-mapping",dest="skip_mapping",default=False,action="store_true",help="Skip mapping")
	group_steps.add_argument("--sorted",dest="sorted_bam",default=False,action="store_true",help="Do not sort bam")
	group_steps.add_argument("--no-index",dest="no_sorted_bam_index",default=False,action="store_true",help="Do not index bam")
	#Files
	group_files = parser.add_argument_group("Files")
	group_files.add_argument("--del-convert",dest="del_convert",default=False,action="store_true",help="Do ot delete the converted fastq files")
	group_files.add_argument("--del-sam",dest="del_sam",default=False,action="store_true",help="Delete the sam files")
	group_files.add_argument("--del-bam",dest="del_bam",default=False,action="store_true",help="Delete the unsorted bam")
	group_files.add_argument("--no-multi-fastq",dest="no_multi_fastq",default=False,action="store_true",help="No multimapping in unmapped fastq")
	#Other
	group_other = parser.add_argument_group("Other")
	group_other.add_argument("--version",action="version",version="%(prog)s 1.0")
	options = parser.parse_args()
	
	signal.signal(signal.SIGINT,signal_handler)
	multimapper_num = 0
	unmapped_num = 0
	unique_num = 0
	A2G_num = 0
	T2C_num = 0
	
	if options.continue_prefix:
		hisat2_prefix = options.continue_prefix
	else:
		hisat2_prefix = "hisat2_" + str(os.getpid())
	fwd_A2G = None
	rev_T2C = None
	
	if options.fwd:
		if options.fwd.endswith(".fastq") == False:
			raise ValueError("Forward file do not end with .fastq")
		fwd_A2G = options.fwd.replace(".fastq",".A2G.fastq")
	if options.rev:
		if options.rev.endswith(".fastq") == False:
			raise ValueError("Reverse file do not end with .fastq")
		rev_T2C = options.rev.replace(".fastq",".T2C.fastq")
	
	#Step1 convert transcripts
	if options.skip_convert == False:
		if options.fwd:
			fwd_read_dict,total_reads_fwd = fastq_converter_1(options.fwd,fwd_A2G,"A2G")
		if options.rev:
			rev_read_dict,total_reads_rev = fastq_converter_1(options.rev,rev_T2C,"T2C")
			if total_reads_fwd != total_reads_rev:
				raise Warning("Read number not equal!")
		else:
			total_reads = total_reads_fwd
	
	#Step2 mapping
	if options.skip_mapping == False:
		sys.stderr.write("[%s]Mapping with hisat2, TEMP prefix: %s\n" % (strftime("%Y-%m-%d %H:%M:%S", time.localtime()),hisat2_prefix))
		hisat2 = options.hisat2_path
		hisat2_parameters = read_parameters(options.hisat2_param)
		hisat2_parameters = mapping_string(hisat2_parameters,fwd_A2G,rev_T2C)
		fwd_index = options.index + "/A2G/" + options.index_prefix + "_A2G"
		rev_index = options.index + "/T2C/" + options.index_prefix + "_T2C"
		
		if options.fwd:
			pool = multiprocessing.Pool(2)
			try:
				pool.apply_async(mapping, args=(fwd_index,"A2G",hisat2,hisat2_parameters,hisat2_prefix,))
				pool.apply_async(mapping, args=(rev_index,"T2C",hisat2,hisat2_parameters,hisat2_prefix,))
				pool.close()
				pool.join()
			finally:
				pool.terminate()
	else:
		if not options.continue_prefix:
			raise Warning("Please provide a SAM prefix.")
	
	if options.skip_convert:
		#Step get C-contained reads
		sys.stderr.write("[%s]Finding non-converted C containing reads...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
		fwd_read_dict,rev_read_dict,total_reads = read_m6A_position(options.fwd,options.rev)
	
	if options.del_convert == True:
		os.remove(fwd_A2G)
		if rev_T2C:
			os.remove(rev_T2C)
	
	#Step Get unique mapping and report unmapped and multimappers
	sys.stderr.write("[%s]Handling SAM outputs...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	if options.fwd and options.rev:
		sam_handler_PE()
	elif options.fwd and not options.rev:
		sam_handler_SE()

	if options.del_sam == True:
		if options.fwd:
			os.remove(hisat2_prefix + ".A2G.sam")
			os.remove(hisat2_prefix + ".T2C.sam")
		# os.remove("hisat2.A2G.sam")
		# os.remove("hisat2.T2C.sam")

	unique_num = A2G_num + T2C_num
	total_reads = unique_num + multimapper_num + unmapped_num
	sys.stderr.write("[%s]Completed successfully:\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	sys.stderr.write(" Total reads: %d\n" % total_reads)
	sys.stderr.write(" Unique mapping: %d (%.3f%%)\n" % (unique_num,100*unique_num/(total_reads+0.0)))
	sys.stderr.write("   A2G: %d (%.2f%%)\n" % (A2G_num,100*A2G_num/(total_reads+0.0)))
	sys.stderr.write("   T2C: %d (%.2f%%)\n" % (T2C_num,100*T2C_num/(total_reads+0.0)))
	sys.stderr.write(" Multiple mapping: %d (%.3f%%)\n" % (multimapper_num,100*multimapper_num/(total_reads+0.0)))
	sys.stderr.write(" Unmapped: %d (%.3f%%)\n" % (unmapped_num,100*unmapped_num/(total_reads+0.0)))
	
	if options.sorted_bam == True:
		sys.stderr.write("[%s]Sorting bam...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
		pysam.sort("-o", options.output+".sorted.bam", options.output+".bam")
		if options.no_sorted_bam_index == False:
			sys.stderr.write("[%s]Indexing bam...\n" % strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
			pysam.index(options.output+".sorted.bam")

	if options.del_bam == True:
		os.remove(options.output+".bam")
