import os, time, random
import multiprocessing
from cython_ext.one_population_ancestral_selection_graph import CythonOnePopulationAncestralSelectionGraph


###############################################################################
# This file is a Python script for simulating genealogies WITHIN a population #
# under various conditions using the Ancestral Selection Graph Model          #
###############################################################################

def aRun(mutation_rate, selection_coefficient, population_size, sample_num):
    while True:
        asg = CythonOnePopulationAncestralSelectionGraph(mutation_rate, selection_coefficient, population_size, sample_num)
        topology = asg.pick_embedded_genealogy()
        genotypes = asg.genotypes
        del asg
        if len(set(genotypes)) >= 2:
            break
    return (topology, genotypes)


def runARepeat(tree_num_in_each_repeat, pool, mutation_rate, selection_coefficient, population_size, sample_num):    
    resultList = []
    for i in range(tree_num_in_each_repeat):
        time.sleep(0.001)
        a = pool.apply_async(aRun, (mutation_rate, selection_coefficient, population_size, sample_num))
        resultList.append(a)
    return resultList


def runAParaList(repeat_num, tree_num_in_each_repeat, paraList, outDir):
    if not os.path.isdir(outDir):
        os.makedirs(outDir)

    resultListList = []
    pool = multiprocessing.Pool(processes=int(0.5*multiprocessing.cpu_count())) # Modify it to adjust the number of processes for running the simulation
    for i in range(repeat_num):
        resultList = runARepeat(tree_num_in_each_repeat, pool, *paraList)
        resultListList.append(resultList)
    pool.close()
    pool.join()

    for i, resultList in enumerate(resultListList):
        outFile = os.path.join(outDir, "repeat_"+str(i+1)+".csv")
        with open(outFile, "w") as outHandle:
            outHandle.write("\t".join(["Topology", "Genotype_1", "Genotype_2", "Genotype_3"])+"\n")
            for result in resultList:
                topology, genotypes = result.get()
                outHandle.write("\t".join([topology.decode()] + [str(g) for g in genotypes])+"\n")


def runAll(resultRoot):
    repeat_num = 100
    tree_num_in_each_repeat = 3000 # Each repeat contain 3000 trees to count the distribution of genealogies
    population_size = 2e5
    mutation_rate = 3e-8
    selection_coefficient_list = [-7.5e-6, -5e-6, -2.5e-6, -1e-6, 0]
    sample_num = 3
    for selection_coefficient in selection_coefficient_list:
        outDir = os.path.join(resultRoot, "s="+str(selection_coefficient))
        paraList = [mutation_rate, selection_coefficient, population_size, sample_num]
        runAParaList(repeat_num, tree_num_in_each_repeat, paraList, outDir)


if __name__ == "__main__":
    resultRoot = os.path.join("..", "..", "Results", "Simulation_results_within_population") # The output directory of simulation results
    runAll(resultRoot)
    print("OK")


