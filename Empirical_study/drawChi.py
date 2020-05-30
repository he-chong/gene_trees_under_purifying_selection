from __future__ import division
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
import numpy as np

def makeRegText(slope, intercept, p_value):
    i = 1
    while p_value * 10**i < 10:
        a = p_value * 10**i
        t = 10**i
        i += 1

    p_str = str(round(a,1))+r"${\times} 10^{-"+str(i-1)+"}$" if i >= 3 else str(round(p_value,2))
    text = '   ='+str(round(slope,2))+r"$\mathit{x}$+"+str(round(intercept,2))+", $\mathit{p}$="+p_str if intercept > 0 else '   ='+str(round(slope,2))+r"$\mathit{x}$"+str(round(intercept,2))+", $\mathit{p}$="+p_str
    return text


def drawCSVChi(csvFile, ax, color):
    df = pd.read_csv(csvFile, delimiter="\t")
    textList = []
    keyList = df.keys()[1:4]

    exp = (df["Total"] - df[keyList[-1]])/2
    chi = ((df[keyList[0]]-exp)**2+(df[keyList[1]]-exp)**2)/exp
    chi.name = "Chi"
    df = pd.concat([df, chi], axis=1)

    df.plot(x="Mean_omega", y="Chi", ax=ax, marker="o", markerfacecolor=color, markeredgecolor=color, linestyle="")
    slope, intercept, r_value, p_value, std_error = stats.linregress(df["Mean_omega"], chi)
    reg_x = np.array(df["Mean_omega"])
    reg_y = slope*reg_x + intercept
    ax.plot(reg_x, reg_y, linestyle="--", color=color)

    text = makeRegText(slope, intercept, p_value)

    x_down, x_up = ax.get_xlim()
    y_down, y_up = ax.get_ylim()
    y_up = y_up+(y_up-y_down)*0.46
    ax.set_ylim([y_down, y_up])
    d = 0.08
    ax.text(x_down * 0.95 + x_up * 0.05, y_down * 0.12 + y_up * 0.88,
            text, fontsize=8, color='black')
    ax.text(x_down * 0.95 + x_up * 0.05, y_down * 0.12 + y_up * 0.88,
            r'$\mathit{y}$', fontsize=8, color=color)
    ax.set_xlabel("")
    ax.legend()


def draw(rootDir, outFig, color, nrows, ncols):
    with PdfPages(outFig) as pdf:
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols,figsize=(30,90), sharex=False, sharey=False)
        tripleList = os.listdir(rootDir)
        textListList = []
        for i, tripleBase in enumerate(tripleList):
            print tripleBase
            rowNow = int(i / ncols)
            colNow = i % ncols
            ax = axes[rowNow][colNow]
            # ax.set_title(tripleBase)
            tripleDir = os.path.join(rootDir, tripleBase)
            csvFile = os.path.join(tripleDir, "stat.csv")
            textList = drawCSVChi(csvFile, ax, color)
            textListList.append(textList)
        
        pdf.savefig()
        plt.subplots_adjust(hspace=0.5)

def drawOverallChi(rootDir, outFig, color, groupList):
    with PdfPages(outFig) as pdf:
        fig, ax = plt.subplots(1, 1, figsize=(8,6))
        dfList = []
        for tripleBase in os.listdir(rootDir):
            print tripleBase
            tripleDir = os.path.join(rootDir, tripleBase)
            csvFile = os.path.join(tripleDir, "stat.csv")
            df = pd.read_csv(csvFile, delimiter="\t")
            df.columns = ['Bin_name'] + groupList + ['Total', "Pseudo_Chi", 'Mean_omega']
            dfList.append(df)
        totalDf = pd.concat(dfList)
        totalDf['Mean_omega'] = 1-totalDf['Mean_omega']
        totalDf.plot(x="Mean_omega", y="Pseudo_Chi", ax=ax, marker=".", markersize=1.2, markerfacecolor=color, markeredgecolor=color, linestyle="", legend=False)
        slope, intercept, r_value, p_value, std_error = stats.linregress(totalDf["Mean_omega"], totalDf["Pseudo_Chi"])
        reg_x = np.arange(min(totalDf["Mean_omega"]), max(totalDf["Mean_omega"]), 0.01)
        reg_y = slope*reg_x + intercept
        ax.plot(reg_x, reg_y, linestyle="-", color=color, linewidth=2.5)

        text = makeRegText(slope, intercept, p_value)

        x_down, x_up = ax.get_xlim()
        y_down, y_up = ax.get_ylim()
        ax.text(x_down * 0.95 + x_up * 0.05, y_down * 0.12 + y_up * 0.88,
                text, fontsize=10, color='black')
        ax.text(x_down * 0.95 + x_up * 0.05, y_down * 0.12 + y_up * 0.88,
                r'$\mathit{y}$', fontsize=10, color=color)

        ax.set_xlabel("")
        pdf.savefig()
        plt.subplots_adjust(hspace=0.5)


if __name__ == '__main__':
    groupList = ["Primate", "Scandentia", "Glires"]
    rootDir = os.path.join("..","..","Results","Empirical_results","_".join(groupList),"IV_Statistics_results")
    overallFig = os.path.join("..","..","Results","Empirical_results","_".join(groupList),"statistics_overall_chi.pdf")
    color = "coral"
    drawOverallChi(rootDir, overallFig, color, groupList)