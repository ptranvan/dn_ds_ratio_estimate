# -*- coding: utf-8 -*-

####################################################################################
### Author: Patrick Tran van (patrick.tranvan@unil.ch)
###
### This script is useful for running PAML's codeml with RAxML's tree.
###
### Need these folders and files :
###
### alignments/ : contain the multiple alignments NON formatted in phylip format.
### branch_lengths/ : contain the RAxML's tree NON formatted in newick format.
### codeml.ctl : control template.
###
### At the end, a .sh file is created for VITAL-IT infrastructure. 
###
### How to use:
###
### python mites2codeml.py -s sai_ti -i1 alignments -i2 branch_lengths -i3 codeml_sai_ti.ctl -e <your_email> -o sai_ti.sh
### python mites2codeml.py -s hydrophobicity -i1 alignments -i2 branch_lengths -i3 codeml_hydrophobicity.ctl -e <your_email> -o hydrophobicity.sh
###
####################################################################################

import os
import re
import argparse
import sys
import glob
import os.path
from Bio import SeqIO
import shutil

def alignments_format(alignment_before, alignment_after):
	"""
	Add spaces after header (required for codeml). 
	"""
	
	alignment_before_file = open(alignment_before,"r")
	content = alignment_before_file.readlines()
	
	alignment_after_file = open(alignment_after,"a")			
	alignment_after_file.write(content[0])

	specie_row = 1
	
	while specie_row < len(content):
		alignment_after_file.write(content[specie_row].split()[0] + "   " + content[specie_row].split()[1] + "\n")
		specie_row+=1
	
	alignment_before_file.close()
	alignment_after_file.close()

def codeml_control(organism, model, sample_name, treefile, outfile, dest):
	"""
	Create a codeml's control file for each job. 
	"""

	if organism == "mite_sai_ti":
	
		replacements = {"alignment.phy":sample_name, "mite_sai_ti.tree":treefile, "paml.out":outfile}
		
		with open("codeml_sai_ti.ctl") as infile, open(dest, 'w') as outfile:
		    for line in infile:
		        for src, target in replacements.iteritems():
		            line = line.replace(src, target)
		        outfile.write(line)

	if organism == "mite_hydrophobicity":
	
		replacements = {"alignment.phy":sample_name, "mites_hydrophobicity.tree":treefile, "paml.out":outfile}
		
		with open("codeml_hydrophobicity.ctl") as infile, open(dest, 'w') as outfile:
		    for line in infile:
		        for src, target in replacements.iteritems():
		            line = line.replace(src, target)
		        outfile.write(line)
		        
		        
def modelSAI(raxml_before, raxml_after):
	"""
	Format to SAI model. 
	"""
	
	SAI_before_file = open(raxml_before,"r")
	
	# Read the newick file
	
	for raxml_tree in SAI_before_file:
		
		split_raxml = raxml_tree.split(",")

		complete_format = []
		
		# The branch length are multiplied by 3 and branch label is added (sex or asex)
		
		part0 = split_raxml[0].split(":")[0] + ":" + str(round(float(split_raxml[0].split(":")[1])*3,9)) + "#2" 
		part1 = split_raxml[1].split(":")[0] + ":" + str(round(float(split_raxml[1].split(":")[1])*3,9)) + "#2"
		part2 = split_raxml[2].split(":")[0] + ":" + str(round(float(split_raxml[2].split(":")[1])*3,9))
		part3 = split_raxml[3].split(":")[0] + ":" + str(round(float(split_raxml[3].split(":")[1])*3,9)) + "#2"
		part4 = split_raxml[4].split(":")[0] + ":" + str(round(float(split_raxml[4].split(":")[1].split(")")[0])*3,9)) + "):" +  str(round(float(split_raxml[4].split(":")[2].split(")")[0])*3,9)) + "#1):" + str(round(float(split_raxml[4].split(":")[3].split(")")[0])*3,9)) + "#1):" + str(round(float(split_raxml[4].split(":")[4].split(")")[0])*3,9)) + "#1"
		part5 = split_raxml[5].split(":")[0] +  ":" + str(round(float(split_raxml[5].split(":")[1].split(")")[0])*3,9)) + ");"
		
		complete_format.append(part0)
		complete_format.append(part1)
		complete_format.append(part2)
		complete_format.append(part3)
		complete_format.append(part4)
		complete_format.append(part5)
		
		new_tree = ','.join(complete_format)
				
		SAI_after_file = open(raxml_after,"w")		
		SAI_after_file.write(new_tree)					
		SAI_after_file.close()
		
	SAI_before_file.close()

def modelTI(raxml_before, raxml_after):
	"""
	Format to TI model. 
	"""
	
	TI_before_file = open(raxml_before,"r")
	
	# Read the newick file
	
	for raxml_tree in TI_before_file:
		
		split_raxml = raxml_tree.split(",")

		complete_format = []
		
		# The branch length are multiplied by 3 and branch label is added (sex or asex)
		
		part0 = split_raxml[0].split(":")[0] + ":" + str(round(float(split_raxml[0].split(":")[1])*3,9)) 
		part1 = split_raxml[1].split(":")[0] + ":" + str(round(float(split_raxml[1].split(":")[1])*3,9))
		part2 = split_raxml[2].split(":")[0] + ":" + str(round(float(split_raxml[2].split(":")[1])*3,9))
		part3 = split_raxml[3].split(":")[0] + ":" + str(round(float(split_raxml[3].split(":")[1])*3,9))
		part4 = split_raxml[4].split(":")[0] + ":" + str(round(float(split_raxml[4].split(":")[1].split(")")[0])*3,9)) + "):" +  str(round(float(split_raxml[4].split(":")[2].split(")")[0])*3,9)) + "#1):" + str(round(float(split_raxml[4].split(":")[3].split(")")[0])*3,9)) + "#1):" + str(round(float(split_raxml[4].split(":")[4].split(")")[0])*3,9)) + "#1"
		part5 = split_raxml[5].split(":")[0] +  ":" + str(round(float(split_raxml[5].split(":")[1].split(")")[0])*3,9)) + ");"

		
		complete_format.append(part0)
		complete_format.append(part1)
		complete_format.append(part2)
		complete_format.append(part3)
		complete_format.append(part4)
		complete_format.append(part5)
		
		new_tree = ','.join(complete_format)
		
		TI_after_file = open(raxml_after,"w")		
		TI_after_file.write(new_tree)					
		TI_after_file.close()
		
	TI_before_file.close()

def modelHydro(raxml_before, raxml_after):
	"""
	Format to Hydrophobicity model. 
	"""
	
	Hydrophobicity_before_file = open(raxml_before,"r")
	
	# Read the newick file
	
	for raxml_tree in Hydrophobicity_before_file:
		
		split_raxml = raxml_tree.split(",")
		
		complete_format = []
		
		# The branch length are multiplied by 3
		
		part0 = split_raxml[0].split(":")[0] + ":" + str(round(float(split_raxml[0].split(":")[1])*3,9)) 
		part1 = split_raxml[1].split(":")[0] + ":" + str(round(float(split_raxml[1].split(":")[1])*3,9))
		part2 = split_raxml[2].split(":")[0] + ":" + str(round(float(split_raxml[2].split(":")[1])*3,9))
		part3 = split_raxml[3].split(":")[0] + ":" + str(round(float(split_raxml[3].split(":")[1])*3,9))
		part4 = split_raxml[4].split(":")[0] + ":" + str(round(float(split_raxml[4].split(":")[1])*3,9))
		part5 = split_raxml[5].split(":")[0] + ":" + str(round(float(split_raxml[5].split(":")[1].split(")")[0])*3,9)) + "):" +  str(round(float(split_raxml[5].split(":")[2].split(")")[0])*3,9)) + "):" + str(round(float(split_raxml[5].split(":")[3].split(")")[0])*3,9)) + "):" + str(round(float(split_raxml[5].split(":")[4].split(")")[0])*3,9)) + "):" + str(round(float(split_raxml[5].split(":")[5].split(")")[0])*3,9))
		part6 = split_raxml[6].split(":")[0] + ":" +  str(round(float(split_raxml[6].split(":")[1].split(")")[0])*3,9)) + "):" + split_raxml[6].split(":")[2]
		
		complete_format.append(part0)
		complete_format.append(part1)
		complete_format.append(part2)
		complete_format.append(part3)
		complete_format.append(part4)
		complete_format.append(part5)
		complete_format.append(part6)

		new_tree = ','.join(complete_format)
				
		Hydrophobicity_after_file = open(raxml_after,"w")		
		Hydrophobicity_after_file.write(new_tree)					
		Hydrophobicity_after_file.close()
		
	Hydrophobicity_before_file.close()
				
def sai_ti(alignment_directory, branch_lengths, control_file, email, output_file):
	"""
	Creating script for the SAI and TI models. 
	"""
	
	folder = "sai_ti_output"
	
	if os.path.isfile(output_file):
			os.remove(output_file)
	
	script = open(output_file,"a")
	script.write("#!/bin/bash\n\n")
	script.close()
	
	# Reset the result folder if existing
					
	if os.path.isdir(folder):
		shutil.rmtree(folder, ignore_errors=True)
	else:
		os.mkdir(folder)
	
	# Look for every alignment file and create all the files needed for codeml
	
	alignments_repo = os.path.join(alignment_directory, "*")
	
	alignment_file = 0
	while alignment_file < len(glob.glob(alignments_repo)):
	
		sample = sorted(glob.glob(alignments_repo))[alignment_file]	# alignments file
		sample_name = sample.split("/")[-1]
				
		sample_repo = os.path.join(folder, sample_name)
		modelSAI_repo = os.path.join(sample_repo, "modelSAI")
		modelTI_repo = os.path.join(sample_repo, "modelTI")
	
		# Create a directory for each sample and model
		
		os.makedirs(sample_repo)
		
		os.makedirs(modelSAI_repo)
		alignments_format(sample, os.path.join(modelSAI_repo, sample_name))	# alignment re-format
	
		codeml_control("mite_sai_ti", "2", sample_name, "RAxML_result." + sample_name, sample_name + "SAI", os.path.join(modelSAI_repo, "codeml.ctl"))	# control creation
		modelSAI(os.path.join(branch_lengths, "RAxML_result." + sample_name), os.path.join(modelSAI_repo, "RAxML_result." + sample_name))
		
		os.makedirs(modelTI_repo)
		alignments_format(sample, os.path.join(modelTI_repo, sample_name))	# alignment re-format
		
		codeml_control("mite_sai_ti", "2", sample_name, "RAxML_result." + sample_name, sample_name + "TI", os.path.join(modelTI_repo, "codeml.ctl"))	# control creation
		modelTI(os.path.join(branch_lengths, "RAxML_result." + sample_name), os.path.join(modelTI_repo, "RAxML_result." + sample_name))
	
	
		# Writing a script for VITAL-IT.
		
		script = open(output_file,"a")
				
		script.write("echo 'cd {} && module add Phylogeny/paml/4.9a && codeml' | bsub -u {} -N -J paml_{}\n\n".format(modelSAI_repo, email, sample_name + "_SAI"))
		script.write("echo 'cd {} && module add Phylogeny/paml/4.9a && codeml' | bsub -u {} -N -J paml_{}\n\n".format(modelTI_repo, email, sample_name + "_TI"))
							
		script.close()				
			
		alignment_file += 1

def hydrophobicity(alignment_directory, branch_lengths, control_file, email, output_file):
	"""
	Creating script for the hydrophobicity model. 
	"""
	
	folder = "hydrophobicity_output"
	
	if os.path.isfile(output_file):
			os.remove(output_file)
	
	script = open(output_file,"a")
	script.write("#!/bin/bash\n\n")
	script.close()
			
	# Reset the result folder if existing
					
	if os.path.isdir(folder):
		shutil.rmtree(folder, ignore_errors=True)
	else:
		os.mkdir(folder)
	
	# Look for every alignment file and create all the files needed for codeml
	
	alignments_repo = os.path.join(alignment_directory, "*")
	
	alignment_file = 0
	while alignment_file < len(glob.glob(alignments_repo)):
	
		sample = sorted(glob.glob(alignments_repo))[alignment_file]	# alignments file
		sample_name = sample.split("/")[-1]
				
		sample_repo = os.path.join(folder, sample_name)
		modelHydro_repo = os.path.join(sample_repo, "modelHydrophobicity")
	
		# Create a directory for each sample and model
		
		os.makedirs(sample_repo)
		
		os.makedirs(modelHydro_repo)
		alignments_format(sample, os.path.join(modelHydro_repo, sample_name))	# alignment re-format
	
		codeml_control("mite_hydrophobicity", "2", sample_name, "RAxML_result." + sample_name, sample_name + "hydro", os.path.join(modelHydro_repo, "codeml.ctl"))	# control creation
		modelHydro(os.path.join(branch_lengths, "RAxML_result." + sample_name), os.path.join(modelHydro_repo, "RAxML_result." + sample_name))
			
		# Writing a script for VITAL-IT.
		
		script = open(output_file,"a")
				
		script.write("echo 'cd {} && module add Phylogeny/paml/4.9a && codeml' | bsub -u {} -N -J paml_{}\n\n".format(modelHydro_repo, email, sample_name + "_Hydro"))
							
		script.close()				
			
		alignment_file += 1

##############
# Main script
						
def main(argv):
	
	mod=[]
	mod.append('\n%(prog)s -s sai_ti -i1 <alignments_directory> -i2 <branch_lengths_directory> -i3 <control_file> -e <email> -o <output_file>')
	mod.append('%(prog)s -s hydrophobicity -i1 <alignments_directory> -i2 <branch_lengths_directory> -i3 <control_file> -e <email> -o <output_file>')
	
	parser = argparse.ArgumentParser(prog = 'mites2codeml.py',
                                 usage = "\n".join(mod))

	parser.add_argument('-s', action='store', dest='step_value',
	                    help='Step')
	                                                 	
	parser.add_argument('-i1', action='store', dest='input_value',
	                    help='Input 1')

	parser.add_argument('-i2', action='store', dest='input2_value',
	                    help='Input 2')

	parser.add_argument('-i3', action='store', dest='input3_value',
	                    help='Input 3')

	parser.add_argument('-e', action='store', dest='email_value',
	                    help='Email')
       	
	parser.add_argument('-o', action='store', dest='output_value',
	                    help='Output')
	                    	                    	
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')

	results = parser.parse_args()
	
	# Run PAML codeml in Vital-IT cluster.
	
	if results.step_value == "sai_ti" and results.input_value and results.input2_value and results.input3_value and results.email_value and results.output_value:
		sai_ti(results.input_value, results.input2_value, results.input3_value, results.email_value, results.output_value)

	if results.step_value == "hydrophobicity" and results.input_value and results.input2_value and results.input3_value and results.email_value and results.output_value:
		hydrophobicity(results.input_value, results.input2_value, results.input3_value, results.email_value, results.output_value)
		
											
if __name__ == "__main__":
	main(sys.argv[1:])
				
