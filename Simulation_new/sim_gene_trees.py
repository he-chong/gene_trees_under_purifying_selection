import os, time, random
import multiprocessing
from cython_ext.three_species_ancestral_selection_graph import CythonThreeSpeciesAncestralSelectionGraph

#######################################################################
# This file is a Python script for simulating three-species trees     #
# under various conditions using the Ancestral Selection Graph Model  #
#######################################################################

def aRun(mutation_rate, selection_coefficient, terminal_population_size_1, terminal_population_size_2, terminal_population_size_3, \
            internal_population_size_23, internal_population_size_1_23, divergence_1, divergence_2):
    while True:
        three_species_asg = CythonThreeSpeciesAncestralSelectionGraph(mutation_rate, selection_coefficient, terminal_population_size_1, terminal_population_size_2, terminal_population_size_3, \
            internal_population_size_23, internal_population_size_1_23, divergence_1, divergence_2)
        topology = three_species_asg.pick_embedded_genealogy()
        genotypes = three_species_asg.genotypes
        del three_species_asg
        if len(set(genotypes)) >= 2:
            break
    return (topology, genotypes)

def runARepeat(tree_num_in_each_repeat, pool, mutation_rate, selection_coefficient, terminal_population_size_1, terminal_population_size_2, terminal_population_size_3, \
            internal_population_size_23, internal_population_size_1_23, divergence_1, divergence_2):
    
    resultList = []
    for i in range(tree_num_in_each_repeat):
        time.sleep(random.randrange(1, 3, 1)/1000.0)
        a = pool.apply_async(aRun, (mutation_rate, selection_coefficient, terminal_population_size_1, terminal_population_size_2, terminal_population_size_3, \
            internal_population_size_23, internal_population_size_1_23, divergence_1, divergence_2))
        resultList.append(a)
    return resultList


def runAParaList(repeat_num, tree_num_in_each_repeat, paraList, outDir):
    if not os.path.isdir(outDir):
        os.makedirs(outDir)

    resultListList = []
    pool = multiprocessing.Pool(processes=int(0.5*multiprocessing.cpu_count()))  # Modify it to adjust the number of processes for running the simulation
    for i in range(repeat_num):
        time.sleep(random.randrange(1, 3, 1)/1000.0)
        resultList = runARepeat(tree_num_in_each_repeat, pool, *paraList)
        resultListList.append(resultList)
    pool.close()
    pool.join()

    for i, resultList in enumerate(resultListList):
        outFile = os.path.join(outDir, "repeat_"+str(i+1)+".tre")
        with open(outFile, "w") as outHandle:
            outHandle.write("\t".join(["Topology", "Genotype_1", "Genotype_2", "Genotype_3"])+"\n")
            for result in resultList:
                topology, genotypes = result.get()
                outHandle.write("\t".join([topology.decode()] + [str(g) for g in genotypes])+"\n")

def runAll(resultRoot):
    repeat_num = 100
    tree_num_in_each_repeat = 10000 # Each repeat contain 10000 trees to count the distribution of genealogies
    mutation_rate = 3e-8
    selection_coefficient_list = [-7.5e-6, -5e-6, -2.5e-6, -1e-6, 0] # Set various selection coefficients to the species tree
    internal_population_sizes = [2e5, 2e5]  # Set various population combinations to the terminal branches of the species tree (1, (2, 3))
    terminal_population_sizes_list = [
        [2e5, 2e5, 2e5],
        [2e5, 2e5, 4e4],
        [1e6, 2e5, 2e5],
        [1e6, 2e5, 4e4], 
        [1e6, 1e6, 2e5], 
        [1e6, 1e6, 4e4],
    ]    
    divergence_time_1 = 1.6e6
    divergence_time_2_list = [14e3, 8e3, 2e3] # Set various delta_tau

    for divergence_time_2 in divergence_time_2_list:
        dir_level1 = os.path.join(resultRoot, "delta_tau="+str(divergence_time_2))
        for selection_coefficient in selection_coefficient_list:
            dir_level2 = os.path.join(dir_level1, "s="+str(selection_coefficient))
            for i, terminal_population_sizes in enumerate(terminal_population_sizes_list):
                paraList = [mutation_rate, selection_coefficient] + terminal_population_sizes + internal_population_sizes + [divergence_time_1, divergence_time_2]
                paraStr = "_".join([str(i+1)] + [str(para) for para in paraList])
                dir_level3 = os.path.join(dir_level2, paraStr)
                if not os.path.isdir(dir_level3):
                    runAParaList(repeat_num, tree_num_in_each_repeat, paraList, dir_level3)
    print("OK")


if __name__ == "__main__":
    resultRoot = os.path.join("..", "..", "Results", "Simulation_results_three_species") # The output directory of simulation results
    runAll(resultRoot)


