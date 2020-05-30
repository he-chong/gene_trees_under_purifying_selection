Please make sure Python 2.7, Cython and C/C++ compiler have been installed. 

src/ contains contains C++ codes for simulating a random three-species gene tree given a species tree and random genealogy within a population
In src/ there are makefiles to build demonstrating executable program. One can build the demostrating programs by the following commands
	cd the/path/of/src/
	make -f Makefile_XXX      (If using MinGW in the Windows systems, one may use mingw32-make rather than make)

To build the extension module for simulating three-species gene trees, run the following command:
    python setup_three_species.py build_ext --inplace

If the module is build successfully, there will be a file generated which is named "sim_three_species_ancestral_selection_graph.so" in Unix-like systems or "sim_three_species_ancestral_selection_graph.pyd" in Windows system.

Similarly, one can use setup_within_population.py to build extension module for simulating the genealogy with a population.

Once the extension module is build successfully, one can modify and run the Python script sim_gene_trees.py to generate gene trees under various conditions.