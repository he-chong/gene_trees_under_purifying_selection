import os, re, collections
from Bio import Phylo
from io import StringIO


# This is for compling the results generated by sim_gene_trees.py analyzing the extant genotypes of the three species.


def compileInfo(paraListDir, outFile):
    print(paraListDir)
    with open(outFile, "w") as outHandle:
        outHandle.write("\t".join(["Repeat", "Topology_1", "Topology_2", "Topology_3"])+"\n")
        for repeat in sorted(os.listdir(paraListDir), key=lambda i:int(i.strip(".tre").split("_")[1])):
            repeatFile = os.path.join(paraListDir, repeat)
            genoList1 = []
            genoList2 = []
            genoList3 = []
            with open(repeatFile) as repeatHandle:
                for line in repeatHandle.readlines()[1:]:
                    line_split = line.strip().split("\t")
                    geno1 = line_split[1]
                    geno2 = line_split[2]
                    geno3 = line_split[3]
                    genoList1.append(geno1)
                    genoList2.append(geno2)
                    genoList3.append(geno3)
            count1 = collections.Counter(genoList1)
            count2 = collections.Counter(genoList2)
            count3 = collections.Counter(genoList3)
            countStrList = [str(each["1"]) for each in [count1, count2, count3]]
            outHandle.write("\t".join([repeat.split("_")[1]]+countStrList)+"\n")

                
def compileAll(rootDir):
    for divergence_time_2 in os.listdir(rootDir):
        dir_level1 = os.path.join(rootDir, divergence_time_2)
        if os.path.isdir(dir_level1):
            for selection_coefficient in os.listdir(dir_level1):
                dir_level2 = os.path.join(dir_level1, selection_coefficient)
                if os.path.isdir(dir_level2):
                    for paraList in os.listdir(dir_level2):
                        paraListDir = os.path.join(dir_level2, paraList)
                        if os.path.isdir(paraListDir):
                            outFile = os.path.join(dir_level2, "geno#"+paraList+".csv")
                            compileInfo(paraListDir, outFile)
                

if __name__ == '__main__':
    resultRoot = os.path.join("..", "..", "..", "Results", "Simulation_results_three_species") # The path of the simulation results generated by sim_gene_trees.py
    compileAll(resultRoot)
    print("OK")

