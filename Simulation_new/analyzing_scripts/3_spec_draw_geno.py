from __future__ import division
import os, glob, copy
import csv
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
import numpy as np
import pandas as pd


# This is for drawing the results generated by sim_gene_trees.py
# Before running this script, run compile_geno.py first


def drawForATop(dfForATop, positions, ax, color):
    width=0.2
    meanList = [np.mean(dfForATop[key]) for key in dfForATop]
    stdList = [np.std(dfForATop[key], ddof=1) for key in dfForATop]
    ax.bar(positions, meanList, width, yerr=stdList, color=color, edgecolor="", ecolor="#202020", capsize=2, error_kw={"elinewidth":1})
    ax.plot(positions, meanList, linestyle="", marker='.', markersize=3, color="#202020")

def draw(divergence_time_2_dir, outFig, axes, colorList, firstline):
    fontsize = 10
    for i, selection_coefficient in enumerate(sorted(os.listdir(divergence_time_2_dir), key=lambda i:abs(float(i.split("=")[-1])))):
        ax = axes[i]
        subdir = os.path.join(divergence_time_2_dir, selection_coefficient)
        rawDfList = []
        sort_lambada = lambda f:os.path.basename(f).split("_")[0]
        for statFile in sorted(glob.glob(os.path.join(subdir, "geno#*.csv")),key=sort_lambada):
            df = pd.read_csv(statFile, delimiter='\t')
            total = pd.Series()
            keyList = df.keys()[1:]
            for key in keyList:
                series = df[key]
                if total.empty:
                    total = copy.deepcopy(series)
                else:
                    total += series
            for key in keyList:
                df[key] = df[key]/10000
            rawDfList.append(df)
        topoList = rawDfList[0].keys()[1:]
        for k, topo in enumerate(topoList):
            seriesList = []
            for raw_df in rawDfList:
                series = raw_df[topo]
                seriesList.append(series)
            dfForATop = pd.concat(seriesList, axis=1, keys=[chr(each+97).upper() for each in range(len(rawDfList))])
            positions = np.arange(len(dfForATop.keys()))+1+(k-1)*0.26
            if k==0:
                if firstline:
                    ax.set_title(str(selection_coefficient), fontsize=fontsize)
                ax.set_xticks(range(1,len(rawDfList)+1))
                ax.set_xticklabels([chr(65+x) for x in range(len(rawDfList))])
                ax.xaxis.set_tick_params(labelsize=fontsize)
                ax.yaxis.set_tick_params(labelsize=fontsize)
                ax.set_xlim(0.5, len(rawDfList)+0.5)
                ax.set_ylim(0, 1)
                # ax.set_ylim(0.28, 0.41)
                for x in range(int(round(len(dfForATop/2)))):
                    ax.axvspan(2*x+0.5, 2*x+1.5, facecolor='#f2f2f2', alpha=1)
            drawForATop(dfForATop, positions, ax, colorList[k])
    # plt.show()


def drawAll(rootDir, colorList):
    outFig = os.path.join(rootDir, "result_geno.pdf")
    divergence_time_2_list = []
    for divergence_time_2_base in sorted(os.listdir(rootDir)):
        divergence_time_2_dir = os.path.join(rootDir, divergence_time_2_base)
        if os.path.isdir(os.path.join(rootDir, divergence_time_2_base)):
            divergence_time_2_list.append(divergence_time_2_dir)
    divergence_time_2_list.sort(key=lambda i:float(os.path.basename(i).split("=")[-1]))
    with PdfPages(outFig) as pdf:        
        fig, axes = plt.subplots(nrows=len(divergence_time_2_list), ncols=5, figsize=(13,7), sharex=True, sharey=True)
        for i, divergence_time_2_dir in enumerate(divergence_time_2_list):
            firstline = i == 0
            draw(divergence_time_2_dir, outFig, axes[i], colorList, firstline)
        plt.subplots_adjust(wspace=0.1)
        pdf.savefig()


if __name__ == '__main__':
    resultRoot = os.path.join("..", "..","..","Results", "Simulation_results_three_species")  # The path of the simulation results generated by sim_gene_trees.py
    colorList = [
        "#75b84f",
        "cornflowerblue",
        "#ffad0f",
    ]
    drawAll(resultRoot, colorList)
    print("OK")