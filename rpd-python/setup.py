from setuptools import setup, find_packages

with open("../README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="rpd",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="Python implementation of Ratio Percentile Deviation (RPD) method",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/rpd",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.6",
    install_requires=[
        "numpy>=1.18.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "jupyter>=1.0.0",
            "matplotlib>=3.0",
            "seaborn>=0.10",
            "pandas>=1.0",
        ],
    },
)