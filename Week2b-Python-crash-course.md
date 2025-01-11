## Introduction to Python

**By the end of this module you should be able to**

1. Explain the different data types and structures used in R Programming
  
2. Create functions to analyze your own datasets
  

**Due to limited class time I have summarized the most important materials from these sources. On your own time, I highly recommend you read through the following tutorials**

- [Software Carpentries - R Novice](Programming with R: Summary and Setup]([Programming with R: Summary and Setup](https://swcarpentry.github.io/r-novice-inflammation/))
  

---

### Lesson Materials

- [Python Inflammation Data](https://swcarpentry.github.io/python-novice-inflammation/data/python-novice-inflammation-data.zip)
- [Python Inflammation Code](https://swcarpentry.github.io/python-novice-inflammation/files/code/python-novice-inflammation-code.zip)

---

## First, lets Download neccessary files we need

Navigate to your desktop and create a folder called `swc-python` . Afterwads, use the `curl` commnad to download the two zip file above in the `swc-python` directory and unzip them using the `unzip` command.

```bash
cd Desktop
mkdir -p eDNA-course/swc-python
cd eDNA-course/swc-python
curl -O "https://swcarpentry.github.io/python-novice-inflammation/data/python-novice-inflammation-data.zip""
curl -O "https://swcarpentry.github.io/python-novice-inflammation/files/code/python-novice-inflammation-code.zip""
unzip *.zip
```

## Getting Started with Conda

- Conda is an environment management system that allows you to
  
  - Identify and install compatible versions of software and the required dependencies
    
  - Handle the process of updating software to the most recent versions
    
- It is availabe on Windows, Mac, and Linux machines
  

## CONDA vs Miniconda vs Anaconda

- **Conda** is the management system to instill software packages
  
- **Miniconda** combines conda with Python and a few core packages
  
- **Anaconda** includes miniconda and a large variety of well maintained python packages
  

## What is a Conda Environment

A [Conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html) is a directory that contains a specific collection of Conda packages that you have installed. For example, you may be working on a research project that requires NumPy 1.18 and its dependencies, while another environment associated with an finished project has NumPy 1.12 (perhaps because version 1.12 was the most current version of NumPy at the time the project finished). If you change one environment, your other environments are not affected. You can easily activate or deactivate environments, which is how you switch between them.

**Conda has a default environment called `base` that include a Python installation and some core system libraries and dependencies of Conda. It is a “best practice” to avoid installing additional packages into your `base` software environment. Additional packages needed for a new project should always be installed into a newly created Conda environment.**

## Installing a new Conda environment

Make a new directory where you are going to install this enviroment

```bash
pwd 
```

```bash
/Users/alejandrodesantiago/Desktop/
```

```bash
mkdir conda-env/jupyter
```

create a conda environment in that directory

```bash
conda create -p /Users/alejandrodesantiago/Desktop/conda-env/jupyter
```

```bash
conda activate /Users/alejandrodesantiago/Desktop/eDNA-course/python.3-13
```

Activate the conda environment using the `full path` and then install your software packages - in this case Jupyter

```bash
conda install jupyter
```

## Initializing a Jupyter Notebook

Navigate to the directory where your data files are located:

```bash
cd /Users/alejandrodesantiago/Desktop/eDNA-course/swc-python
```

Use the following command to initialize a jupyter notebook. This will open a browser window on your local-host.

```bash
jupyter notebook
```
