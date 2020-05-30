from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(
    ext_modules=cythonize(
        Extension(
            'sim_three_species_ancestral_selection_graph',
            sources = ["src/cython_three_species_ancestral_selection_graph.pyx", "src/three_species_ancestral_selection_graph.cpp", "src/terminal_ancestral_selection_graph.cpp", \
            "src/internal_ancestral_selection_graph.cpp", "src/base_ancestral_selection_graph.cpp", "src/for_ancestral_selection_graph.cpp"],
            language="c++",
            extra_compile_args=["-std=c++11"]
            )
        )
    )