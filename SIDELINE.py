#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os
import scanpy as sc
import random
from scipy.stats import pearsonr
import networkx as nx
import time
from datetime import timedelta


def checkGeneList(geneList, normCD, allOutput_df):
    for gene in geneList:
        if gene not in normCD.columns:
            print(gene +" not in matrix")
            geneList.remove(gene)
            allOutput_df.drop(gene, inplace=True, axis=1)
    return geneList, allOutput_df


def normalize(df, zeroThresh):
    df.replace(0,np.nan, inplace=True)
    df=df.dropna(thresh=len(df)*zeroThresh, axis=1)
    norm_df = (df - df.min())/(df.max()-df.min())
    norm_df.replace(np.nan, 0, inplace = True)
    return norm_df



def pseudotimeAssign(norm_df):
    adata = sc.AnnData(norm_df, dtype = np.float64)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution = 1.0)
    sc.tl.paga(adata, groups = 'leiden')
    df_res = adata.obs['leiden'].to_frame()
    return df_res



def rankPseudotime(time_df):
    time_df['Rank'] = time_df['leiden'].rank(method = 'first')
    time_df = time_df.sort_values(by = ['Rank'])
    cell_order = list(time_df.index)
    return cell_order



def reorderCells(df, cellOrderList): 
    df = df.T
    df = df[cellOrderList]
    df = df.T
    return df



def pearson_corr(normCD):
    corrM=normCD.corr(method='pearson')
    corrM=corrM.rename_axis('temp')
    result = corrM.stack().reset_index()
    result.columns = ['fromGene','toGene','corr']
    return result



def predict_GRN(geneList, corrResult, numTopRanked, normCD, cudaUse):
    tcdfInput = "./TCDF_Output/inputTemp.csv"
    tcdfOutput= "./TCDF_Output/outputTemp.csv"
    allInteractions = pd.DataFrame(columns = geneList)
    for j in range(len(geneList)):
        interactionsPP = pd.DataFrame(columns = ('gene.pair', 'timesteps'))
        df = corrResult[corrResult['toGene'] == geneList[j]] 
        absList = list()
        for value in df['corr']:
            absolute = abs(value)
            absList.append(absolute)
        df['AbsVal'] = absList 
        topGenes = df.sort_values(['AbsVal'], ascending=[False])
        gm = topGenes['fromGene'].head(numTopRanked + 1)
        gl = gm.tolist()
        geneMatrix = normCD.filter(gl)
        geneMatrix = geneMatrix.replace(np.nan, 0)
        geneMatrix.to_csv(tcdfInput, index=False)
        os.environ['MKL_THREADING_LAYER'] = 'GNU'
        if cudaUse == True:
            command = 'python ./utils/runTCDF.py --data ' + tcdfInput + ' --cuda > ./TCDF_Output/tcdfRunning.txt'    #If you want to use SIDELINE from a different working directory, change this to the absolute path
        else:
            command = 'python ./utils/runTCDF.py --data ' + tcdfInput + ' > ./TCDF_Output/tcdfRunning.txt'
        os.system(command)
        TCDFRes = pd.read_csv(tcdfOutput)
        interactionsPP=interactionsPP.append(TCDFRes, ignore_index=True)
        interactionsPP=interactionsPP[interactionsPP['timesteps'] ==0]
        interactionsPP=interactionsPP.drop('timesteps', axis = 1)
        if len(interactionsPP) > 0:
            interactionsPP[['fromGene','toGene']] = interactionsPP.iloc[:,0].str.split('>',expand=True)
            interactionsPP.insert(len(interactionsPP.T), "Direction", np.nan)
            for r in range(len(interactionsPP)):
                fromGene = normCD.loc[:, interactionsPP.iloc[r,1]]
                toGene = normCD.loc[:, interactionsPP.iloc[r,2]]
                corr, _ = pearsonr(fromGene, toGene)
                if corr > 0:
                    interactionsPP.iloc[r, 3] = "+"
                else: 
                    interactionsPP.iloc[r, 3] = "-" 
            interactionsPP["gene.pair"] = interactionsPP["gene.pair"].astype(str) +","+ interactionsPP["Direction"]
            interactionsPP.drop(columns=['fromGene', 'toGene', 'Direction'], inplace = True)

        t = 0
        while t < len(interactionsPP):
            item = interactionsPP.iloc[t, interactionsPP.columns.get_loc('gene.pair')]
            allInteractions = allInteractions.append({geneList[j]:item}, ignore_index = True)
            t+=1
        print("Gene:\t" + geneList[j])
    return allInteractions

        
        
        
def GRN_background(geneList, numTopGenes, normCD, cudaUse, outDir):
    tcdfInput = "./TCDF_Output/inputTemp.csv"
    tcdfOutput= "./TCDF_Output/outputTemp.csv"
    allInteractions = pd.DataFrame(columns = geneList)
    for j in range(len(geneList)):
        interactionsPP = pd.DataFrame(columns = ('gene.pair', 'timesteps')) 
        cols = list(normCD.columns)
        randomGenes = random.sample(cols, numTopGenes-1)
        randomGenes.insert(0,geneList[j])
        geneMatrix = normCD.filter(randomGenes)        
        geneMatrix = geneMatrix.replace(np.nan, 0)
        os.environ['MKL_THREADING_LAYER'] = 'GNU'
        os.chdir('..')
        geneMatrix.to_csv(tcdfInput, index=False)
        if cudaUse == True:
            command = 'python ./utils/runTCDF.py --data ' + tcdfInput + ' --cuda > /tmp/tcdfRunning.txt'
        else: 
            command = 'python ./utils/runTCDF.py --data ' + tcdfInput + ' > /tmp/tcdfRunning.txt'
        os.system(command)
        TCDFRes = pd.read_csv(tcdfOutput)
        os.chdir(outDir)
        interactionsPP=interactionsPP.append(TCDFRes, ignore_index=True)
        interactionsPP=interactionsPP[interactionsPP['timesteps'] ==0]
        interactionsPP=interactionsPP.drop('timesteps', axis = 1)
        for t in range(len(interactionsPP)):
            item = interactionsPP.iloc[t, interactionsPP.columns.get_loc('gene.pair')]
            allInteractions = allInteractions.append({geneList[j]:item}, ignore_index = True)
        print("Gene:\t" + geneList[j])
    return allInteractions




def combine_permutations(allInteractions, perms, geneList):
    for m in range(len(geneList)):
        counts = allInteractions[geneList[m]].value_counts()  
        outputFileName = 'to' + geneList[m] + '_' + str(perms) + '_counts.csv' 
        counts.to_csv(outputFileName,index=True)  
        
        os.chdir("./Graphs")
        from matplotlib import pyplot as plt
        df = counts.to_frame()
        plt.hist(np.array(df.iloc[:,0]))
        plt.yscale('log')
        plt.title(df.columns[0])
        outputGraphName = 'to' + geneList[m] + '_' + str(perms) + '_counts.png'
        plt.savefig(outputGraphName, dpi=300, bbox_inches='tight') 

        os.chdir('..')
        
        
        
        
def calSigLvl(bg_df, gene): 
    sig_lvl = bg_df.loc[:, gene].mean()
    return sig_lvl



def removeBackground(actual_df, sig_val,gene): 
    actual_noNoise = actual_df[actual_df[gene] >= sig_val]
    actual_noNoise.columns = ['genePair', 'valueCounts']
    return actual_noNoise



def graphGRN(actual_df, gene): 
    actual_df.columns = ['genePair', 'weight']
    actual_df[['genePair','arrowShape']] = actual_df.iloc[:,0].str.split(',',expand=True)
    actual_df[['source','target']] = actual_df.iloc[:,0].str.split('>',expand=True)
    actual_df = actual_df[['source', 'target', 'weight', 'arrowShape']]
    actual_df.replace("+", "delta", inplace=True)
    actual_df.replace("-", "T", inplace=True)

    G = nx.from_pandas_edgelist(actual_df, edge_attr=True)
    fileout = str(gene) + '.graphml'
    nx.write_graphml(G, fileout)

    return actual_df 







def SIDELINE(outDir, permutations, cellNumber, rawCase, rawCtrl, zeroThresh, geneList, allInteractions, topGenes, cudaUse):
    os.mkdir(outDir)
    for i in range(permutations):
        start = time.time()
        case=rawCase.sample(n=cellNumber,axis='columns')  
        ctrl=rawCtrl.sample(n=cellNumber, axis='columns') 
        case=case.T
        ctrl=ctrl.T
        
        
        # 1. normalization + pseudotime
        case = normalize(case, zeroThresh)
        ctrl = normalize(ctrl, zeroThresh) 
        time_case = pseudotimeAssign(case)
        time_ctrl = pseudotimeAssign(ctrl)
        case_cellOrder = rankPseudotime(time_case)
        ctrl_cellOrder = rankPseudotime(time_ctrl)
        case = reorderCells(case, case_cellOrder)
        ctrl = reorderCells(ctrl, ctrl_cellOrder)
        
        
        # 2. difference matrix
        case.reset_index(drop=True, inplace=True)
        ctrl.reset_index(drop=True, inplace=True)
        caseMctrl = case.subtract(ctrl)
        caseMctrl = normalize(caseMctrl, zeroThresh)
        if len(caseMctrl.columns) >= 10000: 
            print("You have " + str(len(caseMctrl.columns)) + " genes after normalization and filtering. You may want to increase the zero threshold for faster performance.")
        else: 
            print(str(len(caseMctrl.columns)) + " genes after normalization and filtering")
        geneList, allInteractions = checkGeneList(geneList, caseMctrl, allInteractions)


        # 3. correlation 
        corr = pearson_corr(caseMctrl)
        
        
        # 4. Detect gene interactions
        interactions = predict_GRN(geneList, corr, topGenes, caseMctrl, cudaUse)
        allInteractions = allInteractions.append(interactions, ignore_index=True)
        del(interactions)

        exec_time = (time.time()-start)
        print("Permutation: \t" + str(i+1))
        print("Run time:\t" + str(timedelta(seconds=exec_time)))


    #combine permutations  
    os.chdir(outDir)
    os.mkdir("Graphs")
    combine_permutations(allInteractions, permutations, geneList)
    

    
    
    

def BACKGROUND(outDir, geneList, permutations, rawCase, rawCtrl, cellNumber, zeroThresh, topGenes, cudaUse): 
    if os.path.isdir(outDir) == True: 
        os.chdir(outDir)
        os.mkdir('./Background')
    else:
        os.mkdir('./Background')
    
    print("Running Background")
    backgroundInteractions = pd.DataFrame(columns = geneList)

    for k in range(permutations):
        start=time.time()
        case=rawCase.sample(n=cellNumber,axis='columns')  
        ctrl=rawCtrl.sample(n=cellNumber, axis='columns')  
        case=case.T
        ctrl=ctrl.T
        
        
        # 1. normalization + pseudotime
        case = normalize(case, zeroThresh)
        ctrl = normalize(ctrl, zeroThresh) 
        time_case = pseudotimeAssign(case)
        time_ctrl = pseudotimeAssign(ctrl)
        case_cellOrder = rankPseudotime(time_case)
        ctrl_cellOrder = rankPseudotime(time_ctrl)
        case = reorderCells(case, case_cellOrder)
        ctrl = reorderCells(ctrl, ctrl_cellOrder)
        
        
        # 2. difference matrix
        case.reset_index(drop=True, inplace=True)
        ctrl.reset_index(drop=True, inplace=True)
        caseMctrl = case.subtract(ctrl)
        caseMctrl = normalize(caseMctrl, zeroThresh)
        
        
        # Removed correlation step for background    
        
        # 4. Detect gene interactions
        interactions = GRN_background(geneList, topGenes, caseMctrl, cudaUse, outDir)
        backgroundInteractions = backgroundInteractions.append(interactions, ignore_index=True)
        
        exec_time = (time.time()-start)
        print("Permutation: \t" + str(k+1))
        print("Time elapsed:\t" + str(timedelta(seconds=exec_time)))


    #combine permutations  
    os.chdir("./Background")
    os.mkdir("Graphs")
    combine_permutations(backgroundInteractions, permutations, geneList)
    
    
    
    
def GRAPH(outDir, geneList, permutations):
    if os.path.isdir(outDir) == True: 
        os.chdir(outDir)
        os.mkdir('./GRN_Visualization')
        os.chdir('./GRN_Visualization')
    else:
        os.mkdir('../GRN_Visualization')
        os.chdir('../GRN_Visualization')
    completeInteractions = pd.DataFrame(columns = ['source', 'target', 'weight', 'arrowShape'])
    
    for g in geneList: 
        actual = pd.read_csv('../to' + g + '_' + str(permutations) + '_counts.csv')
        background = pd.read_csv('../Background/to' + g + '_' + str(permutations) + '_counts.csv')
        sig_ct = calSigLvl(background, g)
        
        actual = removeBackground(actual, sig_ct,g)
        actual = graphGRN(actual, g)
        completeInteractions = completeInteractions.append(actual, ignore_index = True)
        
    G = nx.from_pandas_edgelist(completeInteractions, edge_attr=True)
    fileout = 'SIDELINE_OverallInteractions.graphml'
    nx.write_graphml(G, fileout)
