from setuptools import setup, find_packages

setup(
    name="rinq",
    version="0.1",
    author="Shah Ishmam Mohtashim",
    author_email="smohtas@ncsu.edu",
    description="Quantum-enhanced centrality tool for protein networks",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "networkx",
        "biopython",
        "matplotlib",
        "scipy",
        "dimod",
        "dwave-ocean-sdk"
    ],
)
