from __future__ import division
import os, glob, subprocess, multiprocessing, re, shutil
from Bio import AlignIO
import csv
import random


def preparePartitionFile(alnDir, codon="123", codonPartition=False):
    alnFile = os.path.join(alnDir, "alignment.fasta")
    aln = AlignIO.read(alnFile, "fasta")
    alnLen = aln.get_alignment_length()
    partitionFile = os.path.join(alnDir, "partition.txt")
    with open(partitionFile, "w") as partitionHandle:
        if codonPartition == False:
            partitionHandle.write("DNA, CDS=1-"+str(alnLen))
        else:
            if codon == "123":
                for i in range(3):
                    partitionHandle.write("DNA, CDS%d=%d-%d\\3\n" % (i+1, i+1, alnLen))
            elif codon == "12":
                for i in range(2):
                    partitionHandle.write("DNA, CDS%d=%d-%d\\2\n" % (i+1, i+1, alnLen))

def estimateGeneTree(alnDir, outDir, codon='123', codonParition=False):
    if not os.path.isdir(alnDir):
        raise IOError("No alignment directory")
    if not os.path.isdir(outDir):
        os.makedirs(outDir)

    alnFile = os.path.join(alnDir, "alignment.fasta")
    partitionFile = os.path.join(alnDir, "partition.txt")
    preparePartitionFile(alnDir, codon=codon, codonPartition=codonParition)
    os.chdir(alnDir)
    p = subprocess.Popen(["raxmlHPC-AVX", "-s", alnFile, "-q",partitionFile, "-w", outDir, "-f", "a", "-x", str(random.randint(1,999999)), "-p", str(random.randint(1,999999)),"-#", "200", "-m", "GTRGAMMA", "-n", "tre"], stdout=subprocess.PIPE)
    p.communicate()

def getOmega(mlcFile):
    with open(mlcFile) as mlcHandle:
        lines = mlcHandle.readlines()

    dNList = []
    dSList = []
    for line in lines:
        dNPattern = re.compile(r"dN\s=\s\d*\.?\d*")
        dSPattern = re.compile(r"dS\s=\s\d*\.?\d*")
        try:
            dN = float(dNPattern.search(line).group().split("=")[-1].strip())            
            dS = float(dSPattern.search(line).group().split("=")[-1].strip())

            dNList.append(dN)
            dSList.append(dS)
        except AttributeError:
            continue
    try:
        # if dNList and dSList:
        omega = (sum(dNList)/len(dNList))/(sum(dSList)/len(dSList))
        return omega
    except:
        # print mlcFile
        # print dNList
        # print dSList
        return None

def getAllOmega(omegaDir):
    if not os.path.isdir(omegaDir):
        raise IOError(omegaDir+" does not exists")

    omegaDict= {}
    for omegaBase in os.listdir(omegaDir):
        mlcFile = os.path.join(omegaDir, omegaBase, "mlc")
        omega = getOmega(mlcFile)
        if omega != None:
            omegaDict.update({omegaBase:omega})
        else:
            # print os.path.join(omegaDir, omegaBase)
            pass
    return omegaDict

def estimateGeneTrees(alnRoot, omegaDir, outRoot, codon='123', codonParition=False):
    if not os.path.isdir(alnRoot):
        raise IOError("No alignment directory: "+alnRoot)
    if not os.path.isdir(omegaDir):
        raise IOError("No omega directory: "+omegaDir)
    if os.path.isdir(outRoot):
        return False

    omegaDict = getAllOmega(omegaDir)
    filteredOmega = []
    for item in omegaDict.items():
        if item[1] < 1:
            filteredOmega.append(item)
        else:
            print item

    pool = multiprocessing.Pool(processes=32)
    for alnBase, omega in filteredOmega:
        alnDir = os.path.join(alnRoot, alnBase)
        outDir = os.path.join(outRoot, alnBase)
        pool.apply_async(estimateGeneTree, (alnDir, outDir, codon, codonParition))
        # estimateGeneTree(alnDir, outDir, codon, codonParition)
        # count += 1
        # if count >= size:
        #   count = 0
        #   num += 1
    pool.close()
    pool.join()


if __name__ == '__main__':
    alnDir = r'E:\HC\Neutral_Phylogeny\I_Alignments\Bat_Horse_Cattle_Hedgehog_123\ENSG00000000003_TSPAN6_filtered_NT_2N'
    estimateGeneTree(alnDir, codon='123', codonParition=True)
