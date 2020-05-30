from __future__ import division
import pandas as pd
import os, re, glob
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages 


# This is for compling and draw the results generated by sim_asg.py

def getOuterInner(treeStr):
    stripped = re.sub('\)$', '', re.sub('^\(', '', treeStr.strip().strip(";")))
    elements = stripped.split(',')
    outer_inner = {
        "outer": None,
        "inner": [],
    }
    for elem in elements:
        if "(" not in elem and ")" not in elem:
            outer_inner["outer"] = str(int(elem.strip().split(':')[0]))
        else:
            outer_inner["inner"].append(re.sub('\)$', '', re.sub('^\(', '', elem.strip().split(':')[0])).strip())
    return outer_inner


def statARepeat_ad(csvFile):
    df = pd.read_csv(csvFile, delimiter="\t")
    inGenoList = []
    outGenoList = []
    for index, row in df.iterrows():
        Topology = row["Topology"]
        outer_inner = getOuterInner(Topology)
        outGenoList.append(str(row["Genotype_"+outer_inner["outer"]]))
        inGenoList.append(str(row["Genotype_"+outer_inner["inner"][0]]))
        inGenoList.append(str(row["Genotype_"+outer_inner["inner"][1]]))
    return outGenoList.count("1"), inGenoList.count("1")/2


def statARepeat_disad(csvFile):
    df = pd.read_csv(csvFile, delimiter="\t")
    inGenoList = []
    outGenoList = []
    for index, row in df.iterrows():
        Topology = row["Topology"]
        outer_inner = getOuterInner(Topology)
        outGenoList.append(str(row["Genotype_"+outer_inner["outer"]]))
        inGenoList.append(str(row["Genotype_"+outer_inner["inner"][0]]))
        inGenoList.append(str(row["Genotype_"+outer_inner["inner"][1]]))
    return outGenoList.count("2"), inGenoList.count("2")/2


def statAParaList_ad(resultDir, outFig, ax):
    outDict = {
        "Wild_type_outer": [],
        "Wild_type_inner": [],
    }
    for csvFile in glob.glob(os.path.join(resultDir, "*.csv")):
        outCount, inCount = statARepeat_ad(csvFile)
        outDict["Wild_type_outer"].append(outCount)
        outDict["Wild_type_inner"].append(inCount)
    df = pd.DataFrame(outDict)
    df.plot(marker="o", markersize=3, ax=ax)
    ax.legend(loc=1)


def statAParaList_disad(resultDir, outFig, ax):
    outDict = {
        "Mutant_outer": [],
        "Mutant_inner": [],
    }
    for csvFile in glob.glob(os.path.join(resultDir, "*.csv")):
        outCount, inCount = statARepeat_disad(csvFile)
        outDict["Mutant_outer"].append(outCount)
        outDict["Mutant_inner"].append(inCount)
    df = pd.DataFrame(outDict)
    df.plot(marker="o", markersize=3, ax=ax)
    ax.legend(loc=1)


def statAll(resultRoot):
    fontsize = 14
    resultList = []
    for resultBase in os.listdir(resultRoot):
        resultDir = os.path.join(resultRoot, resultBase)
        if os.path.isdir(resultDir):
            resultList.append(resultDir)
    resultList.sort(key=lambda i: float(os.path.basename(i).split("=")[-1]), reverse=True)
    outFig_ad = os.path.join(resultRoot, "result_ad.pdf")
    with PdfPages(outFig_ad) as pdf_ad:
        fig, axes = plt.subplots(nrows=len(resultList), ncols=1, sharex=True, sharey=True, figsize=(5, 20))
        for resultDir, ax in zip(resultList, axes):
            print resultDir
            if os.path.isdir(resultDir):
                ax.set_title(os.path.basename(resultDir), fontsize=fontsize)
                ax.xaxis.set_tick_params(labelsize=fontsize)
                ax.yaxis.set_tick_params(labelsize=fontsize)
                statAParaList_ad(resultDir, outFig_ad, ax)

        x_down, x_up = ax.get_xlim()
        y_down, y_up = ax.get_ylim()
        y_up = y_up+(y_up-y_down)*0.23
        ax.set_xlim([x_down, x_up+1])
        ax.set_ylim([y_down, y_up])
        plt.subplots_adjust(top=0.9, bottom=0.1)
        pdf_ad.savefig()

    outFig_disad = os.path.join(resultRoot, "result_disad.pdf")
    with PdfPages(outFig_disad) as pdf_disad:
        fig, axes = plt.subplots(nrows=len(resultList), ncols=1, sharex=True, sharey=True, figsize=(5, 20))
        for resultDir, ax in zip(resultList, axes):
            print resultDir
            if os.path.isdir(resultDir):
                ax.set_title(os.path.basename(resultDir), fontsize=fontsize)
                ax.xaxis.set_tick_params(labelsize=fontsize)
                ax.yaxis.set_tick_params(labelsize=fontsize)
                statAParaList_disad(resultDir, outFig_disad, ax)

        x_down, x_up = ax.get_xlim()
        y_down, y_up = ax.get_ylim()
        y_up = y_up+(y_up-y_down)*0.23
        ax.set_xlim([x_down, x_up+1])
        ax.set_ylim([y_down, y_up])
        plt.subplots_adjust(top=0.9, bottom=0.1)
        pdf_disad.savefig()

if __name__ == '__main__':
    resultRoot = os.path.join("..", "..", "Results", "Simulation_results_within_population") # The path of the simulation results generated by sim_asg.py
    statAll(resultRoot)
    print "OK"