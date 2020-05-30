from __future__ import division
import glob, os
import pandas as pd
import numpy as np
from scipy.stats import stats

def anova(resultGivenPara):
    seriesList = []
    for csvFile in glob.glob(os.path.join(resultGivenPara, "stat#*.csv")):
        df = pd.read_csv(csvFile, delimiter="\t")
        seriesList.append(df["Topology_3"])
    return stats.f_oneway(*seriesList)


def anovaAll(resultRoot):
    with pd.ExcelWriter(os.path.join(resultRoot, "anova.xlsx")) as excelWriter:
        resultListList = []
        delta_tauList = []
        for each in os.listdir(resultRoot):
            if os.path.isdir(os.path.join(resultRoot, each)):
                delta_tauList.append(each)
        for delta_tau in sorted(delta_tauList, key=lambda s:float(s.split("=")[-1])):
            if not os.path.isdir(os.path.join(resultRoot, delta_tau)):
                continue
            resultList = []
            selList = sorted(os.listdir(os.path.join(resultRoot, delta_tau)), key=lambda s:float(s.split("=")[-1]), reverse=True)
            for sel in selList:
                resultGivenPara = os.path.join(resultRoot, delta_tau, sel)
                to_str = lambda stat: "f = "+str(stat[0])+", p = "+str(stat[1]).replace("e", " * 10 ^ ")
                anovaStat = anova(resultGivenPara)
                resultList.append(to_str(anovaStat))
            resultListList.append(resultList)
        statDf = pd.DataFrame(resultListList, columns=map(lambda i:i.split("=")[-1], selList), index=sorted(delta_tauList, key=lambda s:float(s.split("=")[-1])))
        print statDf
        statDf.to_excel(excelWriter)


def detectViolationFromNeutrality(resultGivenPara, num_in_each_repeat):
    to_str = lambda ttest: "t = "+str(ttest[0])+", p = "+str(ttest[1]).replace("e", " * 10 ^ ")
    delta_tau = float(os.path.basename(os.path.dirname(resultGivenPara)).split("=")[-1])
    resultList = []
    for csvFile in glob.glob(os.path.join(resultGivenPara, "stat#*.csv")):
        popSize = float(os.path.basename(csvFile).split("_")[-4])
        f_disconcord = np.exp(-delta_tau/popSize)/3
        f_concord = 1-2*f_disconcord
        df = pd.read_csv(csvFile, delimiter="\t")
        ttest1 = stats.ttest_1samp(df["Topology_1"]/num_in_each_repeat, f_concord)
        ttest2 = stats.ttest_1samp(df["Topology_2"]/num_in_each_repeat, f_disconcord)
        ttest3 = stats.ttest_1samp(df["Topology_3"]/num_in_each_repeat, f_disconcord)        
        result = [to_str(ttest1), to_str(ttest2), to_str(ttest3)] 
        resultList.append(result)
    statDf = pd.DataFrame(resultList, columns=["(g1, (g2, g3))", "(g2, (g1, g3))", "(g3, (g2, g1))"], index=[chr(each+97).upper() for each in xrange(len(resultList))])
    return statDf


def detectAsymmetry(resultGivenPara):
    to_str = lambda ttest: "t = "+str(ttest[0])+", p = "+str(ttest[1]).replace("e", " * 10 ^ ")
    delta_tau = float(os.path.basename(os.path.dirname(resultGivenPara)).split("=")[-1])
    resultList = []
    for csvFile in glob.glob(os.path.join(resultGivenPara, "stat#*.csv")):
        df = pd.read_csv(csvFile, delimiter="\t")
        ttest = stats.ttest_ind(df["Topology_2"], df["Topology_3"])
        resultList.append(to_str(ttest))
    statSeries = pd.Series(resultList, name="t-test result", index=[chr(each+97).upper() for each in xrange(len(resultList))])
    return statSeries


def detectAnomalousTree(resultGivenPara):
    combinationF = glob.glob(os.path.join(resultGivenPara, "stat#6*.csv"))[0]
    df = pd.read_csv(combinationF, delimiter="\t")
    return stats.ttest_ind(df["Topology_1"], df["Topology_3"])


def runAll(resultRoot):
    num_in_each_repeat = 10000
    with pd.ExcelWriter(os.path.join(resultRoot, "ttest_neutral.xlsx")) as excelWriter1, pd.ExcelWriter(os.path.join(resultRoot, "ttest_asym.xlsx")) as excelWriter2:
        delta_tauList = []
        for each in os.listdir(resultRoot):
            if os.path.isdir(os.path.join(resultRoot, each)):
                delta_tauList.append(each)
        for delta_tau in sorted(delta_tauList, key=lambda s:float(s.split("=")[-1])):
            for sel in sorted(os.listdir(os.path.join(resultRoot, delta_tau)), key=lambda s:float(s.split("=")[-1]), reverse=True):
                resultGivenPara = os.path.join(resultRoot, delta_tau, sel)
                statDf = detectViolationFromNeutrality(resultGivenPara, num_in_each_repeat)
                statSeries = detectAsymmetry(resultGivenPara)
                statDf.to_excel(excelWriter1, sheet_name=delta_tau + ", " + sel)
                statSeries.to_excel(excelWriter2, sheet_name=delta_tau + ", " + sel)


def ttestForAnomalousTrees(resultRoot):
    with pd.ExcelWriter(os.path.join(resultRoot, "anomalousF.xlsx")) as excelWriter:
        delta_tau = "delta_tau=8000.0"
        selList = ["-2.5e-06", "-5e-06", "-7.5e-06"]
        resultList = []
        for sel in selList:
            resultGivenPara = os.path.join(resultRoot, delta_tau, "s="+sel)
            stat = detectAnomalousTree(resultGivenPara)
            # to_str = lambda stat: "t = "+str(stat[0])+", p = "+str(stat[1])
            resultList.append([str(stat[0]), str(stat[1]).replace("e", " * 10 ^ ")])
        statDf = pd.DataFrame(resultList, index=selList, columns=["t", "p"])
        statDf.to_excel(excelWriter)
        print statDf




if __name__ == '__main__':
    resultRoot = r"E:\HC\Neutral_Phylogeny\Results\Simulation_results_three_species"
    # anovaAll(resultRoot)
    # runAll(resultRoot)
    ttestForAnomalousTrees(resultRoot)
    print "OK"