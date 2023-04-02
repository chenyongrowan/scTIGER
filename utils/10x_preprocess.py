#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import gzip
import os
import scipy.io
import pandas as pd
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-d", "--directoryPath", dest="dir", help="Directory path for files to convert from 10x to a gene expression matrix",type=str)
args = parser.parse_args()

mat = scipy.io.mmread(os.path.join(args.dir, "matrix.mtx.gz"))
features_path = os.path.join(args.dir, "features.tsv.gz")
feature_ids = [row[0] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
gene_names = [row[1] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
feature_types = [row[2] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
barcodes_path = os.path.join(args.dir, "barcodes.tsv.gz")
barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, mode="rt"), delimiter="\t")]
matrix = pd.DataFrame.sparse.from_spmatrix(mat)
matrix.columns = barcodes
matrix.insert(loc=0, column="gene", value=gene_names)
filename = args.dir + "/geneExpressionMatrix.csv"
matrix.to_csv(filename, index=False)
print('Finished converting the 10x files to a gene expression matrix. File saved to:\t' + filename)