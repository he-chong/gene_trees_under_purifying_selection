This is a new version of the simulation.
In this version, the implement details are modified so that the simulation can run faster and use less memmory.
Meanwhile, the codes of the new version are easier for understanding.
The result of the new version is the same as before.

Please make sure Python 3.7, Cython and C/C++ compiler have been installed. 

The src/ directory contains C++ codes for simulating a random three-species gene tree 
given a species tree and random genealogy within a population.

To build the extension module for simulating three-species gene trees, run the following command:
    python setup.py build_ext --inplace

If the module is build successfully, there will be two files generated in the /cython_ext directory which are 
suffixed with ".so" in Unix-like systems or with ".pyd" in Windows system.

Once the extension module is built successfully, you can modify and run the Python scripts 
sim_gene_trees.py and sim_one_population.py to perform the simulation.

In addition, in src/ directory there are two makefiles for building demonstration programs. 
These demostration programs can be built by the following commands
	cd the/path/of/src/
	make -f Makefile_XXX      (If using MinGW in the Windows systems, one may use mingw32-make rather than make)