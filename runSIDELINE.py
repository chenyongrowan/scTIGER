#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Madison Dautle
02.02.2023
SIDELINE - The SIngle-cell, DEep-Learning gene regulatory network INferenceE method
"""

import os
import argparse
import sys 
import warnings 
import pandas as pd
import SIDELINE_def as sd
import shutil

warnings.simplefilter('ignore')
os.environ['KMP_WARNINGS'] = 'off'

parser = argparse.ArgumentParser()
#Flags 
parser.add_argument("-p", "--permutations", dest = "runs", default = 100, help="Number of permutations to run. Default 100", type=int)
parser.add_argument("-top", '--numTopGenes', dest = 'topGenes', default = 100, help = "Number of top correlated genes selected. Default 100", type=int)
parser.add_argument("-zero", "--zeroThresh", dest = "zeroThresh", default = 0.40, help="Threshold for number of 0's tolerated for a gene. Default 0.40", type=float)
parser.add_argument("-s", "--start", dest = "start", default = 1, help="Starting point for SIDELINE. Default 1 (Run SIDELINE and background). 2 is background only", type=int)
parser.add_argument("-goi", "--geneOfInterest", dest = "goi", help="One or more genes of interest separated by +")
parser.add_argument("-ctrl", "--control", dest = "ctrl", help="Control condition. A csv file with cells as columns and genes as rows. Must contain at least 10 cells.")
parser.add_argument("-exp", "--experimental",dest ="exp", help="Case/experimental condition. A csv file with cells as columns and genes as rows. Must contain at least 10 cells.")
parser.add_argument("--cuda", dest = 'cuda', action = 'store_true', help="CUDA use on. Default is off.")
parser.add_argument("-o", "--output", dest = 'outputDir', default = 'SIDELINE_Output', help ="Output directory name. Default is 'SIDELINE_Output", type=str)
args = parser.parse_args()

#Arg checking 
#Check that files exist then read in files
fileCheck = os.path.isfile(args.exp)
if fileCheck == False:
    sys.exit("Experimental data file does not exist")
caseW = pd.read_csv(args.exp, index_col = 0)
if len(caseW) < 10: 
    sys.exit("Too few cells in the experimental data file. Must be at least 10 cells.")
fileCheck = os.path.isfile(args.ctrl)
if fileCheck == False:
    sys.exit("Control data file does not exist")
ctrlW = pd.read_csv(args.ctrl,index_col = 0)
if len(ctrlW) < 10: 
    sys.exit("Too few cells in the control data file. Must be at least 10 cells.")

#Separate goi
geneList = args.goi.split("+")
print("Genes of interest: " + str(geneList))
if len(geneList) < 1:
    sys.exit("Need to specify at least 1 gene of interest. Separate multiple genes of interest with +. Make sure ")

#zeroThresh
if args.zeroThresh < 0 or args.zeroThresh > 1: 
    sys.exit("Zero threshold needs to be decimal between 0 and 1") 
    
allInteractions = pd.DataFrame(columns = geneList)

#random select the same number of cells for case and control 
if len(caseW.T) < len(ctrlW.T):
    n = len(caseW.T)
else:
    n = len(ctrlW.T)


os.mkdir('./TCDF_Output')
    
#######################################################################################################################################
while args.start <= 3:
    if args.start == 1: 
        sd.SIDELINE(args.outputDir, args.runs, n, caseW, ctrlW, args.zeroThresh, geneList, allInteractions, args.topGenes, args.cuda)
        os.system("touch commandDetails.txt")
        with open("commandDetails.txt", 'a') as file:
           l1 = "# permutations: " + args.runs
           l2 = "\n# top genes: " + args.topGenes
           l3 = "\nZero Threshold: " + args.zeroThresh
           l4 = "\nControl file name: " + args.ctrl
           l5 = "\nExperimental file name: " + args.case
           l6 = "\nCuda use was on: " + args.cuda 
        file.writelines([l1, l2, l3, l4, l5, l6])
        file.close()
        args.start+=1
    elif args.start == 2:
        sd.BACKGROUND(args.outputDir, geneList, args.runs, caseW, ctrlW, n, args.zeroThresh, args.topGenes, args.cuda)
        args.start+=1
    elif args.start == 3:
        sd.GRAPH(args.outputDir, geneList, args.runs)
    else:
        args.start+=1


os.system('pwd')
shutil.rmtree('../../TCDF_Output')
print('SIDELINE finished.')
