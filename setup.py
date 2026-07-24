from setuptools import setup, find_packages

setup(
    name="mtcojo_postgwas",
    version="1.0.0",
    description="End-to-end GCTA mtCOJO, PostGWAS harmonisation, and LDSC reporting pipeline",
    author="JJOHN41",
    packages=find_packages(),
    package_data={"mtcojo_postgwas": ["reporting/assets/*.R"]},
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
