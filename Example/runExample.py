#!/usr/bin/env python3

import os

os.system('./run_scTIGER.py -goi AR+PTEN+ERG -ctrl ./Example/Data/Patient4_Benign_endothelial.csv -exp ./Example/Data/Patient4_Tumor_endothelial.csv -p 50 -top 100 -zero 0.15 -o ./Example/SampleResult_ProstateCancer')
