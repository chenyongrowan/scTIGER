#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os
import scanpy as sc
from scipy.stats import pearsonr
import networkx as nx
import time
from datetime import timedelta
import bambi as bmb
import arviz as az
from scipy import stats 
import math


def checkGeneList(geneList, normCD, allOutput_df):
    for gene in geneList:
        if gene not in normCD.columns:
            print(gene +" not in matrix")
            geneList.remove(gene)
            allOutput_df = allOutput_df.drop(gene, axis=1)
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



def predict_GRN(geneList, corrResult, numTopRanked, normCD, cudaUse,timeDelay):
    tcdfInput = "./TCDF_Output/inputTemp.csv"
    tcdfOutput= "./TCDF_Output/outputTemp.csv"
    allInteractions = pd.DataFrame(columns = geneList)
    for j in range(len(geneList)):
        interactionsPP = pd.DataFrame(columns = ('gene.pair', 'timesteps'))
        df = corrResult[corrResult['toGene'] == geneList[j]] 
        absList = pd.DataFrame({'AbsVal': df['corr'].abs()})
        df = pd.concat([df, absList], axis=1) 
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
        interactionsPP = pd.concat([interactionsPP, TCDFRes], ignore_index=True)
        if len(interactionsPP) > 0:
            interactionsPP[['fromGene','toGene']] = interactionsPP.iloc[:,0].str.split('>',expand=True)
            interactionsPP.insert(len(interactionsPP.T), "Direction", np.nan)
            for r in range(len(interactionsPP)):
                fromGene = normCD.loc[:, interactionsPP.iloc[r,2]]
                toGene = normCD.loc[:, interactionsPP.iloc[r,3]]
                corr, _ = pearsonr(fromGene, toGene)
                if corr > 0:
                    interactionsPP.iloc[r, 4] = "+"
                else: 
                    interactionsPP.iloc[r, 4] = "-" 
            interactionsPP["gene.pair"] = interactionsPP["gene.pair"].astype(str) +","+ interactionsPP["Direction"]
            interactionsPP.drop(columns=['fromGene', 'toGene', 'Direction'], inplace = True)
        if timeDelay ==0:
            interactionsPP=interactionsPP[interactionsPP['timesteps'] ==0]
            interactionsPP=interactionsPP.drop('timesteps', axis = 1)
        else:
            interactionsPP=interactionsPP[interactionsPP['timesteps'] <= timeDelay]
            interactionsPP['gene.pair'] = interactionsPP['gene.pair'] + ',' + interactionsPP['timesteps'].astype(str)
            interactionsPP=interactionsPP.drop('timesteps', axis = 1)

        t = 0
        while t < len(interactionsPP):
            item = interactionsPP.iloc[t, interactionsPP.columns.get_loc('gene.pair')]
            new_df = pd.DataFrame({geneList[j]: [item]})
            allInteractions = pd.concat([allInteractions, new_df], ignore_index=True)
            t+=1
        print("Gene:\t" + geneList[j])
    return allInteractions

        
        


def combine_permutations(allInteractions, perms, geneList):
    for m in range(len(geneList)):
        counts = allInteractions[geneList[m]].value_counts()  
        outputFileName = geneList[m] + '.csv' 
        counts.to_csv(outputFileName,index=True)  
        
        os.chdir("./Graphs")
        from matplotlib import pyplot as plt
        df = counts.to_frame()
        plt.hist(np.array(df.iloc[:,0]))
        plt.yscale('log')
        plt.title(df.columns[0])
        outputGraphName = geneList[m] + '.png'
        plt.savefig(outputGraphName, dpi=300, bbox_inches='tight') 
        plt.clf()
        
        os.chdir('..')
        
        
        
        
def calSigLvl(bg_df, gene, alpha_val): 
    bgAppCts = bg_df[gene].value_counts().sort_index()
    bgAppCts = bgAppCts.to_frame()
    bgAppCts['Freq'] = bgAppCts.index
    try: 
        eq1 = gene+" ~ Freq"
        model = bmb.Model(eq1, bgAppCts, family = "negativebinomial")
        negB = model.fit()
        summary = az.summary(negB)
        z_val = stats.norm.ppf(1-alpha_val)    
        sig_lvl = summary.loc['Intercept', 'mean'] + ((z_val*summary.loc['Intercept', 'sd'])/math.sqrt(len(bg_df)))
    except:
        print("Negative binomial fitting failed. Taking average to calculate significance value.")
        sig_lvl = bg_df.loc[:, gene].mean()
    return sig_lvl


def removeBackground(actual_df, sig_val, gene): 
    actual_noNoise = actual_df[actual_df[gene] >= sig_val]
    return actual_noNoise



def graphGRN(actual_df, gene, timeDelay): 
    actual_df.columns = ['genePair', 'weight']
    try:
        if timeDelay == 0:
            actual_df[['genePair','arrowShape']] = actual_df.iloc[:,0].str.split(',',expand=True)
        else:
            actual_df[['genePair','arrowShape','timeDelay']] = actual_df.iloc[:,0].str.split(',',expand=True)
        actual_df[['source','target']] = actual_df.iloc[:,0].str.split('>',expand=True)
        if timeDelay ==0:
            actual_df = actual_df[['source', 'target', 'weight', 'arrowShape']]
        else:
            actual_df = actual_df[['source', 'target', 'weight', 'arrowShape', 'timeDelay']]
        actual_df.replace("+", "delta", inplace=True)
        actual_df.replace("-", "T", inplace=True)
    
        G = nx.from_pandas_edgelist(actual_df, edge_attr=True)
        fileout = str(gene) + '.graphml'
        nx.write_graphml(G, fileout)
            
            
    except:
        if timeDelay ==0:
            actual_df = pd.DataFrame(columns = ['source', 'target', 'weight', 'arrowShape'])
        else:
            actual_df = pd.DataFrame(columns = ['source', 'target', 'weight', 'arrowShape', 'timeDelay'])
    return actual_df 







def scTIGER(outDir, permutations, cellNumber, rawCase, rawCtrl, zeroThresh, geneList, allInteractions, topGenes, cudaUse, timeDelay):
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
        if len(caseMctrl.columns) >= 10000: 
            print("You have " + str(len(caseMctrl.columns)) + " genes after normalization and filtering. You may want to increase the zero threshold for faster performance.")
        else: 
            print(str(len(caseMctrl.columns)) + " genes after normalization and filtering")
        geneList, allInteractions = checkGeneList(geneList, caseMctrl, allInteractions)


        # 3. correlation 
        corr = pearson_corr(caseMctrl)
        
        
        # 4. Detect gene interactions
        interactions = predict_GRN(geneList, corr, topGenes, caseMctrl, cudaUse, timeDelay)
        allInteractions = pd.concat([allInteractions, interactions], ignore_index=True)
        del(interactions)

        exec_time = (time.time()-start)
        print("Permutation: \t" + str(i+1))
        print("Run time:\t" + str(timedelta(seconds=exec_time)))


    #combine permutations  
    os.chdir(outDir)
    os.mkdir("Graphs")
    combine_permutations(allInteractions, permutations, geneList)
    
    hold = 1
    return hold
 
    
 
    
def GRAPH(outDir, geneList, permutations, timeDelay, alpha_val, hold):
    if hold == 1: 
        os.mkdir('./GRN_Visualization')
        os.chdir('./GRN_Visualization')
    else:
        os.chdir(outDir)
        os.mkdir('./GRN_Visualization')
        os.chdir('./GRN_Visualization')
        
    if timeDelay ==0: 
        completeInteractions = pd.DataFrame(columns = ['source', 'target', 'weight', 'arrowShape'])
    else:
        completeInteractions = pd.DataFrame(columns = ['source', 'target', 'weight', 'arrowShape', 'timeDelay'])

    for g in geneList: 
        actual = pd.read_csv('../' + g + '.csv')
        actual.columns = ['Interactions', g]
        sig_ct = calSigLvl(actual, g, alpha_val)
        
        actual = removeBackground(actual, sig_ct, g)
        actual = graphGRN(actual, g, timeDelay)
        completeInteractions = pd.concat([completeInteractions, actual], ignore_index=True)
        
    G = nx.from_pandas_edgelist(completeInteractions, edge_attr=True)
    fileout = 'OverallInteractions.graphml'
    nx.write_graphml(G, fileout)

