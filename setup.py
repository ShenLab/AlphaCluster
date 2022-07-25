
import sys
import os
import io
from setuptools import setup
from distutils.core import Extension
from Cython.Build import cythonize

EXTRA_COMPILE_ARGS = ["-std=c++2a"]
EXTRA_LINK_ARGS = []

if sys.platform == "darwin":
    EXTRA_COMPILE_ARGS += ["-stdlib=libc++", "-mmacosx-version-min=10.9"]
    EXTRA_LINK_ARGS += ["-stdlib=libc++", "-mmacosx-version-min=10.9"]

weights = cythonize([
    Extension("alphacluster.weights",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        sources=["alphacluster/weights.pyx",
            "src/weighted_choice.cpp",
            "src/simulate.cpp"],
        include_dirs=["src/"],
        language="c++"),
    Extension("alphacluster.transcript",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        sources=[
            "alphacluster/transcript.pyx",
            "src/tx.cpp"],
        include_dirs=["src/"],
        language="c++"),
    Extension("alphacluster.gencode",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        sources=[
            "alphacluster/gencode.pyx",
            "src/gencode.cpp",
            "src/gtf.cpp",
            "src/tx.cpp",
            "src/gzstream/gzstream.C",
        ],
        include_dirs=["src/", "src/gzstream"],
        libraries=["z"],
        language="c++"),    
    Extension("alphacluster.site_specific_rates",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        sources=[
            "alphacluster/site_specific_rates.pyx",
            "src/weighted_choice.cpp",
            "src/tx.cpp",
            "src/site_rates.cpp"],
        include_dirs=["src/"],
        language="c++"),
    ])

setup(name="alphacluster",
    description='Package to examine de novo clustering',
    long_description=io.open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    version="1.0.0",
    author="Joseph Obiajulu, Jeremy McRae",
    author_email="joseph.obiajulu@gmail.com, jeremy.mcrae@gmail.com",
    license="MIT",
    url='https://github.com/ShenLab/MVP3D',
    packages=["alphacluster"],
    install_requires=[
        'aiohttp >= 3.0',
        'scipy >= 0.9.0',
        'cython >= 0.19.0',
        'pyfaidx >= 0.5.8',
    ],
    package_data={"alphacluster": ['data/rates.txt', 'weights.pxd']},
    entry_points={'console_scripts': ['alphacluster = alphacluster.__main__:main']},
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires='>=3.6',
    ext_modules=weights,
    test_suite="tests")
