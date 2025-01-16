## Introduction to Python

**By the end of this module you should be able to**

1. Explain the different data types and structures used in R Programming
  
2. Create functions to analyze your own datasets
  

**Due to limited class time I have summarized the most important materials from these sources. On your own time, I highly recommend you read through the following tutorials**

- [Software Carpentries - Python]([Programming with Python: Summary and Setup](https://swcarpentry.github.io/python-novice-inflammation/)
  
- [GACRC - Using CONDA](https://wiki.gacrc.uga.edu/images/b/b0/Using_conda_on_the_GACRC_Sapelo2_Cluster.pdf)
  

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

### Introduction to Python data types

There are several kinds of data types but the most common are:

1. Integers
  
2. Floating point numbers
  
3. Strings
  

We can assign a value to a variable using the `=` sign. For example, lets track the weight of a patient by assigning the weight in kg, as a floating point, to a variable called `weight_kg`

```python
weight_kg = 60.3
```

We are also going to create a personal identifier for each patient

```python
patient_id = '001'
```

No that we have our data saved as a variable, we can manipulate them and use them in calculations. Let's convert the weight to pounds and add a prefix to our patient identifier.

```python
weight_lb = 2.2 * weight_kg
patient_id = 'inflam_' + patient_id
```

We can use the function `print` to display things on our screen. We are going to use the `print` function to print the patient ID and their weight to our console:

```python
print(patient_id, 'weight in kilograms:', weight_kg
```

```python
inflam_001 weight in kilograms: 60.3
```

### Loading Data into Python

In order to process data in python we need to use a python package called **numpy**. Numpy stands for numerical python and allows us to manipulate matrices and arrays. We can import packages into python using the `import` function:

```python
import numpy
```

now we can import our datafile and using a numpy function and saving it into a variable called `data`

```python
data = numpy.loadtxt(fname='inflammation-01.csv', delimiter=',')
```

```python
print(data)
```

```
[[ 0.  0.  1. ...,  3.  0.  0.]
 [ 0.  1.  2. ...,  1.  0.  1.]
 [ 0.  1.  1. ...,  2.  1.  1.]
 ...,
 [ 0.  1.  1. ...,  1.  1.  1.]
 [ 0.  0.  0. ...,  0.  2.  0.]
 [ 0.  0.  1. ...,  1.  1.  0.]]
```

You can access a single number from the array using an index in square brackets, but be careful. Programming languages like Fortran, MATLAB and R start counting at 1 because that’s what human beings have done for thousands of years. Languages in the C family (including C++, Java, Perl, and Python) count from 0 because it represents an offset from the first value in the array (the second value is offset by one index from the first value). This is closer to the way that computers represent arrays (if you are interested in the historical reasons behind counting indices from zero, you can read [Mike Hoye’s blog post](https://exple.tive.org/blarg/2013/10/22/citation-needed/)). As a result, if we have an M×N array in Python, its indices go from 0 to M-1 on the first axis and 0 to N-1 on the second. It takes a bit of getting used to, but one way to remember the rule is that the index is how many steps we have to take from the start to get the item we want.

Therefore, to access the first numer in the array we need to use the index [0,0]:

```python
print('first value in data:', data[0, 0])
```

```python
first value in data: 0.0
```

### Analyzing data

The `numpy` package has several useful functions to summarize our datasets. We are going to find the max, min, and standard deviation of our datasets and save them in individual variables using a single line

```python
maxval, minval, stdval = numpy.amax(data), numpy.amin(data), numpy.std(data)

print('maximum inflammation:', maxval)
print('minimum inflammation:', minval)
print('standard deviation:', stdval)
```

```python
maximum inflammation: 20.0
minimum inflammation: 0.0
standard deviation: 4.61383319712
```

### Summarizing data from single patient

We can summarize data pertaining to a single patient by slicing our dataset and saving it as a new temporary array:

```python
patient_0 = data[0, :] # 0 on the first axis (rows), everything on the second (columns)
print('maximum inflammation for patient 0:', numpy.amax(patient_0))
```

```python
maximum inflammation for patient 0: 18.0
```

### Visualizing data

We are going to import the python library `matplotlib.pyplot` to plot our data.

```python
import matplotlib.pyplot
```

Now we will plot the average inflammation over time:

```python
ave_inflammation = numpy.mean(data, axis=0)
ave_plot = matplotlib.pyplot.plot(ave_inflammation)
matplotlib.pyplot.show()
```

We can also plot the min and max inflammation using the following commands:

```python
max_plot = matplotlib.pyplot.plot(numpy.amax(data, axis=0))
matplotlib.pyplot.show()

min_plot = matplotlib.pyplot.plot(numpy.amin(data, axis=0))
matplotlib.pyplot.show()
```

When we are publishing our manuscript, we might need to make a figure with multiple subplots. In this case we can create multiple plots to have a publishable multi-panel figure

 The function `matplotlib.pyplot.figure()` creates a space into which we will place all of our plots. The parameter `figsize` tells Python how big to make this space. Each subplot is placed into the figure using its `add_subplot` [method](https://swcarpentry.github.io/python-novice-inflammation/reference.html#method). The `add_subplot` method takes 3 parameters. The first denotes how many total rows of subplots there are, the second parameter refers to the total number of subplot columns, and the final parameter denotes which subplot your variable is referencing (left-to-right, top-to-bottom). Each subplot is stored in a different variable (`axes1`, `axes2`, `axes3`). Once a subplot is created, the axes can be titled using the `set_xlabel()` command (or `set_ylabel()`).

```python
import numpy
import matplotlib.pyplot

data = numpy.loadtxt(fname='inflammation-01.csv', delimiter=',')

fig = matplotlib.pyplot.figure(figsize=(10.0, 3.0))

axes1 = fig.add_subplot(1, 3, 1)
axes2 = fig.add_subplot(1, 3, 2)
axes3 = fig.add_subplot(1, 3, 3)

axes1.set_ylabel('average')
axes1.plot(numpy.mean(data, axis=0))

axes2.set_ylabel('max')
axes2.plot(numpy.amax(data, axis=0))

axes3.set_ylabel('min')
axes3.plot(numpy.amin(data, axis=0))

fig.tight_layout()

matplotlib.pyplot.savefig('inflammation.png')
matplotlib.pyplot.show()
```

### Repeating functions using Loops
```
import glob
import numpy
import matplotlib.pyplot

filenames = sorted(glob.glob('inflammation*.csv'))
filenames = filenames[0:3]
for filename in filenames:
    print(filename)

    data = numpy.loadtxt(fname=filename, delimiter=',')

    fig = matplotlib.pyplot.figure(figsize=(10.0, 3.0))

    axes1 = fig.add_subplot(1, 3, 1)
    axes2 = fig.add_subplot(1, 3, 2)
    axes3 = fig.add_subplot(1, 3, 3)

    axes1.set_ylabel('average')
    axes1.plot(numpy.mean(data, axis=0))

    axes2.set_ylabel('max')
    axes2.plot(numpy.amax(data, axis=0))

    axes3.set_ylabel('min')
    axes3.plot(numpy.amin(data, axis=0))

    fig.tight_layout()
    matplotlib.pyplot.show()
```
