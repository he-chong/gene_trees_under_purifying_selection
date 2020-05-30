from __future__ import division
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
import numpy as np


def makeRegText(slope, intercept, p_value):
    i = 1
    while p_value * 10**i < 10:
        a = p_value * 10**i
        t = 10**i
        i += 1

    if i >= 3:            
        p_str = str(round(a,1))+r"${\times} 10^{-"+str(i-1)+"}$"
    else:
        p_str = str(round(p_value,2))
    text = '    ='+str(round(slope,2))+r"$\mathit{x}$+"+str(round(intercept,2))+", $\mathit{p}$="+p_str
    return text


def drawCSV(csvFile, ax, colorList, order):
    # print csvFile
    df = pd.read_csv(csvFile, delimiter="\t")
    textList = []
    df['Mean_omega'] = 1-df['Mean_omega']
    keyList = df.keys()[1:4]

    for i, color in zip(order, colorList):
        key = keyList[i]
        # print key
        df[key] = df[key]/df["Total"]
        df.plot(x='Mean_omega', y=key, ax=ax, marker="o", markersize=1, markerfacecolor=color, markeredgecolor=color, linestyle="", legend=False)
        slope, intercept, r_value, p_value, std_error = stats.linregress(df["Mean_omega"], df[key])
        reg_x = np.arange(min(df["Mean_omega"]), max(df["Mean_omega"]), 0.01)
        reg_y = slope*reg_x + intercept
        ax.plot(reg_x, reg_y, linestyle="-", color=color, linewidth=0.5)
        text = makeRegText(slope, intercept, p_value)        
        textList.append(text)

    x_down, x_up = ax.get_xlim()
    y_down, y_up = ax.get_ylim()
    y_up = y_up+(y_up-y_down)*0.46
    ax.set_ylim([y_down, y_up])
    x_down, x_up = ax.get_xlim()
    y_down, y_up = ax.get_ylim()
    d = 0.1
    x_rel_loc = 0.05
    y_rel_loc = 0.88
    for k, color in enumerate(colorList):
        ax.text(x_down * (1-x_rel_loc) + x_up * x_rel_loc, y_down * (1-y_rel_loc+d*k) + y_up * (y_rel_loc-d*k),
                textList[k], fontsize=3, color='black')
        ax.text(x_down * (1-x_rel_loc) + x_up * x_rel_loc, y_down * (1-y_rel_loc+d*k) + y_up * (y_rel_loc-d*k),
                r'$\mathit{y_'+str(k+1)+'}$', fontsize=3, color=color)

    ax.set_xlabel("")
    ax.legend(loc=1, fontsize=3)
    ax.xaxis.set_tick_params(labelsize=4)
    ax.yaxis.set_tick_params(labelsize=4)

def draw(rootDir, outFig, colorList, order, nrows, ncols):
    with PdfPages(outFig) as pdf:
        tripleList = os.listdir(rootDir)
        for page in range(int(len(tripleList)/(nrows*ncols))):
            tripletNow = tripleList[page*nrows*ncols:(page+1)*nrows*ncols]
            fig = plt.figure(figsize=(8,15))
            gs = GridSpec(nrows+1, ncols)
            for i, tripleBase in enumerate(tripletNow):
                print tripleBase
                rowNow = int(i / ncols)
                colNow = i % ncols
                ax = fig.add_subplot(gs[rowNow, colNow])
                ax.set_title(tripleBase.replace(".", ", "), fontsize=4)
                tripleDir = os.path.join(rootDir, tripleBase)
                csvFile = os.path.join(tripleDir, "stat.csv")
                drawCSV(csvFile, ax, colorList, order)
            ax = fig.add_subplot(gs[nrows, :])
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['right'].set_visible(False)
            if page == int(len(tripleList)/(nrows*ncols))-1:
                ax.text(0.07, 0.5, "Figure S1. The results of 288 runs. Each run corresponds to a four-taxon combination. Each of the three possible \ntopologies is represented by the taxon in the relative basal position.",\
                {"family":"Times New Roman", "weight":"normal", "size":10}, horizontalalignment='left', verticalalignment='center')
            else:
                ax.text(0.46, 0.5, "Continued",\
                {"family":"Times New Roman", "weight":"normal", "size":10}, horizontalalignment='left', verticalalignment='center')

            plt.subplots_adjust(top=0.9, bottom=0.1, left=0.05, right=0.95, hspace=0.48, wspace=0.4)
            pdf.savefig()
            

def drawOverall(rootDir, outFig, groupList, colorList, markerList, order):
    with PdfPages(outFig) as pdf:
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        ax.set_ylim([0.15, 0.6])
        dfList = []
        for tripleBase in os.listdir(rootDir):
            # print tripleBase            
            # ax.set_title(tripleBase)
            tripleDir = os.path.join(rootDir, tripleBase)
            csvFile = os.path.join(tripleDir, "stat.csv")
            df = pd.read_csv(csvFile, delimiter="\t")
            df.columns = ['Bin_name']+groupList+['Total','Pseudo_Chi', 'Mean_omega']
            dfList.append(df)
        totalDf = pd.concat(dfList)
        print sum(totalDf["Total"])

        textList = []
        totalDf['Mean_omega'] = 1-totalDf['Mean_omega']
        keyList = totalDf.keys()[1:4]
        x = 0
        for i, color in zip(order, colorList):
            key = keyList[i]
            print key
            totalDf[key] = totalDf[key]/totalDf["Total"]
            print totalDf[key]
            totalDf.plot(x='Mean_omega', y=key, ax=ax, marker=markerList[x], markersize=1.5, markerfacecolor=color, markeredgecolor=color, linestyle="", legend=False)
            x += 1

        for i, color in zip(order, colorList):
            key = keyList[i]
            slope, intercept, r_value, p_value, std_error = stats.linregress(totalDf["Mean_omega"], totalDf[key])    
            reg_x = np.arange(min(totalDf["Mean_omega"]), max(totalDf["Mean_omega"]), 0.01)
            reg_y = slope*reg_x + intercept
            ax.plot(reg_x, reg_y, linestyle="-", color=color, linewidth=2)
            text = makeRegText(slope, intercept, p_value)
            textList.append(text)

        x_down, x_up = ax.get_xlim()
        y_down, y_up = ax.get_ylim()
        d = 0.08
        x_rel_loc = 0.05
        y_rel_loc = 0.9
        for j, color in enumerate(colorList):
            ax.text(x_down * (1-x_rel_loc) + x_up * x_rel_loc, y_down * (1-y_rel_loc+d*j) + y_up * (y_rel_loc-d*j),
                    textList[j], fontsize=12, color='black')
            ax.text(x_down * (1-x_rel_loc) + x_up * x_rel_loc, y_down * (1-y_rel_loc+d*j) + y_up * (y_rel_loc-d*j),
                    r'$\mathit{y_'+str(j+1)+'}$', fontsize=12, color=color)
        ax.set_xlabel("")
        ax.legend(loc=1)
        pdf.savefig()



if __name__ == '__main__':
    groupList = ["Primate", "Scandentia", "Glires"]
    rootDir = os.path.join("..","..","Results","Empirical_results","_".join(groupList),"IV_Statistics_results_50")
    outFig = os.path.join("..","..","Results","Empirical_results","_".join(groupList),"statistics_50.pdf")
    overallFig = os.path.join("..","..","Results","Empirical_results","_".join(groupList),"statistics_overall_50.pdf")
    colorList = [
        "#75b84f",
        "#ffad0f",
        "cornflowerblue",
    ]
    markerList = [".", ".", "."]
    order = [2, 0, 1]
    draw(rootDir, outFig, colorList, order, 9, 4)
    # drawOverall(rootDir, overallFig, groupList, colorList, markerList, order)
    