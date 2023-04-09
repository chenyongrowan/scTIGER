#!/usr/bin/env python3

import os

os.chdir(..)
os.system(./runSIDELINE.py -goi AR+PTEN+ERG -ctrl ./Data/ProstateCancer/Patient4_Benign_endothelial.csv -exp ./Data/ProstateCancer/Patient4_Tumor_endothelial.csv -p 50 -top 100 -zero 0.15 -o SIDELINE_SampleResult_ProstateCancer)
