#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Madison Dautle
Created: 02.02.2023
Updated: 05.12.2023
scTIGER - Single-cell Temporal Inference of Gene Regulatory Networks
"""

import os
import argparse
import sys 
import warnings 
import pandas as pd
import scTIGER as sd
import shutil

warnings.simplefilter('ignore')
os.environ['KMP_WARNINGS'] = 'off'

parser = argparse.ArgumentParser()
#Flags 
parser.add_argument("-p", "--permutations", dest = "runs", default = 100, help="Number of permutations to run. Default 100", type=int)
parser.add_argument("-top", '--numTopGenes', dest = 'topGenes', default = 100, help = "Number of top correlated genes selected. Default 100", type=int)
parser.add_argument("-zero", "--zeroThresh", dest = "zeroThresh", default = 0.30, help="Threshold for number of 0's tolerated for a gene. Default 0.40", type=float)
parser.add_argument("-t", "--timesteps", dest = "timeDelay", default = 0, help = "The number of steps allowed between detected interactions. Default is 0", type = int)
parser.add_argument("-s", "--start", dest = "start", default = 1, help="Starting point for scTIGER. Default 1 (Run scTIGER and make GRN files). 2 is only for generating GRN files from existing scTIGER output (if you want to change the alpha value, you can just run step 2).", type=int)
parser.add_argument("-goi", "--geneOfInterest", dest = "goi", help="One or more genes of interest separated by +")
parser.add_argument("-ctrl", "--control", dest = "ctrl", help="Control condition. A csv file with cells as columns and genes as rows. Must contain at least 10 cells.")
parser.add_argument("-exp", "--experimental",dest ="exp", help="Case/experimental condition. A csv file with cells as columns and genes as rows. Must contain at least 10 cells.")
parser.add_argument("--cuda", dest = 'cuda', action = 'store_true', help="CUDA use on. Default is off.")
parser.add_argument("-o", "--output", dest = 'outputDir', default = 'scTIGER_Output', help ="Output directory name. Default is 'scTIGER_Output", type=str)
parser.add_argument("-a", "--alpha", dest = 'alpha', default = 0.05, help='Alpha value to determine significiant interactions threshold. Default 0.05', type=float)
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

#timesteps 
if args.timeDelay > 2:
    print("We don't suggest using more than 2 timesteps between interactions. You might want to decrease the number of timesteps you are allow.")
    
#Alpha
if args.alpha < 0 or args.alpha > 1:
    sys.exit("That alpha value is invalid. It must be between 0 and 1, but 0.05 is recommended.")

allInteractions = pd.DataFrame(columns = geneList)

#random select the same number of cells for case and control 
if len(caseW.T) < len(ctrlW.T):
    n = len(caseW.T)
else:
    n = len(ctrlW.T)

hold = 0
os.mkdir('./TCDF_Output')
    
#######################################################################################################################################
while args.start <= 2:
    if args.start == 1: 
        hold = sd.scTIGER(args.outputDir, args.runs, n, caseW, ctrlW, args.zeroThresh, geneList, allInteractions, args.topGenes, args.cuda, args.timeDelay)
        os.system("touch commandDetails.txt")
        with open('commandDetails.txt', 'a') as file:
           l1 = "# permutations: " + str(args.runs)
           l2 = "\n# top genes: " + str(args.topGenes)
           l3 = "\nZero Threshold: " + str(args.zeroThresh)
           l4 = "\nControl file name: " + args.ctrl
           l5 = "\nExperimental file name: " + args.exp
           l6 = "\nCuda use was on: " + str(args.cuda)
           l7 = "\nThe maximum number of timesteps was: " + str(args.timeDelay)
           l8 = "\nAlpha value was: " + str(args.alpha)
           file.writelines([l1, l2, l3, l4, l5, l6, l7, l8])
        args.start+=1
    elif args.start == 2:
        sd.GRAPH(args.outputDir, geneList, args.runs, args.timeDelay, args.alpha, hold)
        args.start+=1


shutil.rmtree('../../TCDF_Output')
print('scTIGER finished.')
