from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(
	ext_modules=cythonize(
		Extension(
			'sim_ancestral_selection_graph_within_a_population',
			sources = ["src/cython_ancestral_selection_graph.pyx", "src/ancestral_selection_graph.cpp", 
			"src/base_ancestral_selection_graph.cpp", "src/for_ancestral_selection_graph.cpp"],
			language="c++",
			extra_compile_args=["-std=c++11"]
			)
		)
	)