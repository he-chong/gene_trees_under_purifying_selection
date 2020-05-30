from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(
    ext_package="cython_ext",
    ext_modules=cythonize([
        Extension(
            'three_species_ancestral_selection_graph',
            sources = ["pyx/three_species_ancestral_selection_graph.pyx", "src/three_species_ancestral_selection_graph.cpp", "src/ancestral_selection_graph_composer.cpp", \
            "src/ancestral_selection_graph_component.cpp", "src/for_ancestral_selection_graph.cpp"],
            language="c++",
            extra_compile_args=["-std=c++11"]
            ),
        Extension(
            'one_population_ancestral_selection_graph',
            sources = ["pyx/one_population_ancestral_selection_graph.pyx", "src/one_population_ancestral_selection_graph.cpp", "src/ancestral_selection_graph_composer.cpp", \
            "src/ancestral_selection_graph_component.cpp", "src/for_ancestral_selection_graph.cpp"],
            language="c++",
            extra_compile_args=["-std=c++11"]
            )
        ])
    )