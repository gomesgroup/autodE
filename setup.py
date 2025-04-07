from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import platform

if platform.system() == "Windows":  # compile with Visual C/C++
    cpp_compile_args = ["-permissive-", "-O2"]
    cpp_link_args = []
else:  # on Linux or MacOS
    cpp_compile_args = ["-std=c++11", "-Wno-missing-braces", "-O3"]
    cpp_link_args = ["-std=c++11"]

extensions = [
    Extension("cconf_gen", ["autode/conformers/cconf_gen.pyx"]),
    Extension(
        "ade_dihedrals",
        sources=["autode/ext/ade_dihedrals.pyx"],
        include_dirs=["autode/ext/include"],
        language="c++",
        extra_compile_args=cpp_compile_args,
        extra_link_args=cpp_link_args,
    ),
    Extension(
        "ade_rb_opt",
        sources=["autode/ext/ade_rb_opt.pyx"],
        include_dirs=["autode/ext/include"],
        language="c++",
        extra_compile_args=cpp_compile_args,
        extra_link_args=cpp_link_args,
    ),
]

setup(
    ext_modules=cythonize(extensions, language_level="3"),
    package_data={
        "autode.transition_states": ["lib/*.txt"],
        "autode.solvent": ["lib/*.xyz"],
    },
)
