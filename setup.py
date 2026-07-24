from setuptools import setup, find_packages

setup(
    name="mtcojo_postgwas",
    version="1.0.0",
    description="Unified GCTA mtCOJO & PostGWAS Harmonisation Pipeline Package",
    author="JJOHN41",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "polars>=0.20.0",
        "numpy<2",
        "matplotlib",
        "matplotlib-venn",
    ],
    entry_points={
        "console_scripts": [
            "mtcojo-postgwas = mtcojo_postgwas.cli:main",
        ],
    },
)
