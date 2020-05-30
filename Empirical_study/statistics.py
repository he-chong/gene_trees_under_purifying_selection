from __future__ import division
import os, glob, csv, re, multiprocessing
from collections import OrderedDict, Counter
import numpy as np
import sample
from estimateGeneTrees import getAllOmega

SUPPORT_CUTOFF = 70

class TaxonomyInfo(OrderedDict):   
    def __init__(self, infoStrTuple):
        super(TaxonomyInfo, self).__init__(infoStrTuple)
        self.outgroup = infoStrTuple[-1][1]

def stat(geneTreeDir, omegaDir, outDir, taxonomyInfo, concord, size=300):
    class Tree(object):
        def __init__(self, treeFile, omega, taxonomyInfo):
            with open(treeFile) as treeHandle:
                treeStr = treeHandle.readline()
            searched = re.search('(?<=\))\d+', treeStr)
            if isinstance(searched.group(), str):
                support = int(searched.group())
            else:
                raise Exception() 

            processed = re.sub('\)$', '', re.sub('^\(', '', re.sub('\d*\:\d*\.?\d*', '', treeStr.strip().strip(";"))))
            outtaxon = taxonomyInfo.outgroup
            elements = processed.split(',')
            clustered = []
            non_clustered = []    
            for elem in elements:
                if "(" in elem:
                    clustered.append(elem.strip("("))
                elif ")" in elem:
                    clustered.append(elem.strip(")"))
                else:
                    non_clustered.append(elem)
            if outtaxon in non_clustered:
                i = non_clustered.index(outtaxon)
                outer = non_clustered[i+1-len(non_clustered)]
            else:
                i = clustered.index(outtaxon)
                outer = clustered[i+1-len(clustered)]

            self.treeFile = treeFile
            self.treeStr = treeStr
            self.name = os.path.basename(os.path.dirname(treeFile))
            self.outer = outer
            self.support = support
            self.omega = omega

    class Bin(object):
        def __init__(self, treeObjList, name, taxonomyInfo):
            self.trees = treeObjList
            self.name = name
            self.counter = OrderedDict()
            rawCounter = Counter([tree.outer for tree in self.trees])
            for key in taxonomyInfo.keys()[:3]:
                taxon = taxonomyInfo[key]
                self.counter.update({taxon:rawCounter[taxon]})

        @property
        def meanOmega(self):
            return np.mean([treeObj.omega for treeObj in self.trees])


        def writeLoci(self, lociFile):
            with open(lociFile, "w") as lociHandle:
                writer = csv.writer(lociHandle, delimiter="\t")
                for treeObj in self.trees:
                    treeBase = os.path.basename(os.path.dirname(treeObj.treeFile))
                    writer.writerow([treeObj.name, str(treeObj.omega)])

    def filteredByOmega(omegaDir):
        omegaDict = getAllOmega(omegaDir)
        filteredOmega = []
        for item in omegaDict.items():
            if item[1] < 1:
                filteredOmega.append(item)
        sortedOmega = sorted(filteredOmega, key=lambda i: i[1], reverse=True)
        return sortedOmega

#######
    if not os.path.isdir(outDir):
        os.makedirs(outDir)

    sortedOmega = filteredByOmega(omegaDir)
#binning
    binList = []
    lociCount = 0
    binCount = 0
    treeObjList = []
    for treeBase, omega in sortedOmega:
        treeFile = os.path.join(geneTreeDir, treeBase, "RAxML_bipartitions.tre")
        treeObj = Tree(treeFile, omega, taxonomyInfo)
        if treeObj.support > SUPPORT_CUTOFF:
            treeObjList.append(treeObj)
            lociCount += 1
        else:
            continue
        if lociCount >= size:
            aBin = Bin(treeObjList, "Bin_" + str(binCount+1), taxonomyInfo)
            lociFile = os.path.join(outDir, aBin.name + ".csv")
            aBin.writeLoci(lociFile)
            binList.append(aBin)
            lociCount = 0
            binCount += 1
            treeObjList = []
    aBin = Bin(treeObjList, "Bin_" + str(binCount+1), taxonomyInfo)
    lociFile = os.path.join(outDir, aBin.name + ".csv")
    aBin.writeLoci(lociFile)
    binList.append(aBin)
#output
    outFile = os.path.join(outDir, "stat.csv")
    with open(outFile, "w") as outHandle:
        writer = csv.writer(outHandle, delimiter="\t")
        keyList = binList[0].counter.keys()
        writer.writerow(["Bin_name"]+keyList+["Total", "Pseudo_Chi", "Mean_omega"])
        for aBin in binList: 
            valueList = [aBin.counter[key] for key in keyList]
            total = sum(valueList)
            if concord == 0:
                nc1, nc2 = 1, 2
            elif concord == 1:
                nc1, nc2 = 0, 2
            elif concord == 2:
                nc1, nc2 = 0, 1
            else:
                raise Exception("Unknown concord")

            try:
                exp = (total-valueList[-1])/2
                chi = ((valueList[0]-exp)**2+(valueList[1]-exp)**2)/exp         
                writer.writerow([aBin.name] + valueList + [total, chi, aBin.meanOmega])
            except:
                pass


def statAll(sampleInfo, resultRoot, concord):
    keyList = sampleInfo.keys()
    pool = multiprocessing.Pool(processes=5)
    for taxonStr1 in sampleInfo[keyList[0]]:
        for taxonStr2 in sampleInfo[keyList[1]]:
            for taxonStr3 in sampleInfo[keyList[2]]:
                for outtaxonStr in sampleInfo[keyList[3]]:
                    taxaList = [each[0] for each in [taxonStr1, taxonStr2, taxonStr3, outtaxonStr]]
                    tripleStr = ".".join(taxaList)
                    omegaDir = os.path.join(resultRoot, "II_Omega", tripleStr)  
                    geneTreeDir = os.path.join(resultRoot, "III_Gene_Trees", tripleStr)
                    outDir = os.path.join(resultRoot, "IV_Statistics_results", tripleStr)
                    infoStrTuple = zip(keyList, taxaList)
                    taxonomyInfo = TaxonomyInfo(infoStrTuple)
                    # stat(geneTreeDir, omegaDir, outDir, taxonomyInfo, concord)
                    pool.apply_async(stat, (geneTreeDir, omegaDir, outDir, taxonomyInfo, concord))
    pool.close()
    pool.join()


if __name__ == '__main__':
    resultRoot = os.path.join("..", "..", "Results", "Empirical_results", "Primate_Scandentia_Glires")
    statAll(sample.psgInfo, resultRoot, 2)
