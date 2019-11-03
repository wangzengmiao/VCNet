################################################################################################################
# File: Readme.txt
# Aim: A brief introduction about the usage of VCNet.
# Author:  Zengmiao WANG
# Email : wangzengmiao@gmail.com
# Date  : 2016-08-16
#---------------------------------------------------------------------------------------------------------------
# Package required: CompQuadForm,survey
# function: PerformVCNet
# Input: 
#	Datafile------the name of the data file with ".txt" extension.
#	Resultfile----the names of the result file with ".txt" extension.
# Output:
#	A ".txt" file with the name of "Resultfile".
#---------------------------------------------------------------------------------------------------------------
# Usage:
#	In R environment, type the following commands:
#	source("VCNet.R");
#	PerformVCNet("SimulationData.txt","test_result.txt");
# "SimulationData.txt" is the name of data file and "test_result.txt" is the name of result file.
#---------------------------------------------------------------------------------------------------------------	
# The format of data file and the output file:
#	Both of the two files are text file.
#	Input-------An example of data file (SimulationData.txt) is given in the package. In this file, the first 
#				column is the gene ID, the second column is the exon ID, and the expression of samples are 
#				from the third column to the last column. 
#	Output------An example of result file (test_result.txt) is given in the package. In this file, each row 
#				represents an possible edge. The first and second column are the two genes' IDs. The third
#				column is the test statistics of the possible edge and the last column is the p-value.
#---------------------------------------------------------------------------------------------------------------
