## Introduction to the RStudios

**By the end of this module you should be able to**

1. Explain the different data types and structures used in R Programming
  
2. Create functions to analyze your own datasets
  

**Due to limited class time I have summarized the most important materials from these sources. On your own time, I highly recommend you read through the following tutorials**

- [Software Carpentries - R Novice](Programming with R: Summary and Setup](https://swcarpentry.github.io/r-novice-inflammation/)
  
- [Software Carpentries - R for Reproducible Scientific Analysis](R for Reproducible Scientific Analysis: Summary and Setup](https://swcarpentry.github.io/r-novice-gapminder/)
  
- [Advanced R]([4 Subsetting | Advanced R](https://adv-r.hadley.nz/subsetting.html))
  

---

Install R and Rstudio on your computer:

- [Download and install the latest version of R](https://www.r-project.org/).
- [Download and install RStudio](https://www.rstudio.com/products/rstudio/download/#download).

---

## Workflow within RStudio

There are two main ways one can work within RStudio:

1. Test and play within the interactive R console then copy code into a .R file to run later.
  
2. Start writing in a .R file and use RStudio’s short cut keys for the Run command to push the current line, selected lines or modified lines to the interactive R console.
  

## R Packages

It is possible to add functions to R by writing a package, or by obtaining a package written by someone else. As of this writing, there are over 10,000 packages available on CRAN (the comprehensive R archive network). R and RStudio have functionality for managing packages:

- You can see what packages are installed by typing `installed.packages()`
- You can install packages by typing `install.packages("packagename")`, where `packagename` is the package name, in quotes.
- You can update installed packages by typing `update.packages()`
- You can remove a package with `remove.packages("packagename")`
- You can make a package available for use with `library(packagename)`

Packages can also be viewed, loaded, and detached in the Packages tab of the lower right panel in RStudio. Clicking on this tab will display all of the installed packages with a checkbox next to them. If the box next to a package name is checked, the package is loaded and if it is empty, the package is not loaded. Click an empty box to load that package and click a checked box to detach that package.

Packages can be installed and updated from the Package tab with the Install and Update buttons at the top of the tab.

There are several R packages that are commonly used by microbial ecologists and statisticians including:

- Phyloseq: Explore microbiome profiles using R: https://joey711.github.io/phyloseq/
  
- Vegan: https://cran.r-project.org/web/packages/vegan/vegan.pdf
  
- Metacoder: https://grunwaldlab.github.io/metacoder_documentation/
  
- GGPLOT for data visualization: https://ggplot2.tidyverse.org
  
- ggpubr for publication read graphics and stats: https://rpkgs.datanovia.com/ggpubr/
  

We will go through some of these packages throughout the class. However, the goal of today is get everyone up and running on RStudio.

###

### Data types and Structures

In R there are 6 datatypes:

- character
  
- numeric (real or decimal)
  
- integer
  
- logical
  
- complex
  
- raw (bytes - binary data)
  

You can create data structures by combining data types:

- Vector
  
- List
  
- Matrix
  
- Data Frame
  
- Factors
  

Vectors and List require the data structure to have the same datatypes. While lists and dataframes can head multiple data types. A **vector** is a collection of elements. A **matrix** is a vector with dimensions (i.e., columns and rows). A **dataframe** is similar to a matrix (dimensional) but with difference data types. A **list** is the most flexible - it can be describes as a bin with different data types or structures. **Factors** are specifically categorical variables.

We will go through the different datatypes and structures in this module.

Saving Variables

We can store values in variables using the assignment operator `<-` or `=` like this:

```r
x <- 1/40
```

```r
x = 1/40 
```

The right hand side of the assignment can be any valid R expression. The right hand side is *fully evaluated* before the assignment occurs.

Variable names can contain letters, numbers, underscores and periods but no spaces. They must start with a letter or a period followed by a letter (they cannot start with a number nor an underscore). Variables beginning with a period are hidden variables. Different people use different conventions for long variable names, these include

- periods.between.words
- underscores_between_words
- camelCaseToSeparateWords

What you use is up to you, but **be consistent**.

You can also reassign the variable using the same command

```r
x <- 100
```

One final thing to be aware of is that R is *vectorized*, meaning that variables and functions can have vectors as values. In contrast to physics and mathematics, a vector in R describes a set of values in a certain order of the same data type. For example:

```r
x <- 1:5
2^x
```

### Creating a function

We can also create function to make a group of commands into one single function. For example, lets make a function to convert C to F.

```r
fahrenheit_to_celsius <- function(temp_F) {
  temp_C <- (temp_F - 32) * 5 / 9
  return(temp_C)
}
```

We define `fahrenheit_to_celsius` by assigning it to the output of `function`. The list of argument names are contained within parentheses. Next, the [body](https://swcarpentry.github.io/r-novice-inflammation/reference.html#function-body) of the function–the statements that are executed when it runs–is contained within curly braces (`{}`). The statements in the body are indented by two spaces, which makes the code easier to read but does not affect how the code operates.

![loadingag13663](file:///Users/alejandrodesantiago/Downloads/function-terminology.svg?msec=1735354233040?msec=1735699252927)

Notice the structure of a function:

- **Function Name**: In this case the function is saved as `fahrenheit_to_celsius`. We can call this function using `fahrenheit_to_celsius()`
  
- **Formals**: This is the list of arguments that you can control when you call the function. In this case the only argument we can change when we call the function is the `temp_F`. This variable can be any numerical value.
  
- **Body**: This is the code within the brackets `{}`
  
- **Environment**: The data structure that finds the values assocaited each variable or names. Unless you specify a different environment, this will be the R Global Enviornment which is a list of objects created when we start an R Session.
  

Then you can run your function like this:

```r
# freezing point of water
fahrenheit_to_celsius(32)
```

```
[1] 0
```

We can create a new functio nto convert Celsius to Kelvin

```r
celsius_to_kelvin <- function(temp_C) {
  temp_K <- temp_C + 273.15
  return(temp_K)
}

# freezing point of water in Kelvin
celsius_to_kelvin(0)
```

```r
[1] 273.15
```

Now if we wanted to convert Fehrenheit to Kelvin, we can compose multiple functions together

```r
fahrenheit_to_kelvin <- function(temp_F) {
  temp_C <- fahrenheit_to_celsius(temp_F)
  temp_K <- celsius_to_kelvin(temp_C)
  return(temp_K)
}

# freezing point of water in Kelvin
fahrenheit_to_kelvin(32.0)
```

```r
[1] 273.15
```

Alternatively, we can also nest two functions together.

```r
# freezing point of water in Fahrenheit
celsius_to_kelvin(fahrenheit_to_celsius(32.0))
```

### Creating and manipulating dataframe

We can read a data to a variable using the command read.csv.

```r
dat <- read.csv(file = "data/inflammation-01.csv", header = FALSE)
```

The functions `read.csv`has a single mandatory argument - `file`. This is the path to the CSV file that you are trying to read. The header argument is optional but in this case we are telling it that the file we are reading does NOT have column headers. We can display the first few rows using the function `head`

```r
head(data)
```

We can check the type of data structure the vairable is by using the function `class`

```r
class(dat)
```

```r
[1] "data.frame"
```

Futhermore we can see the dimenstion using the function `dim`

```r
dim(dat)
```

```r
[1] 60 40
```

This tells us that our dataframe has 60 rows and 40 columns. If e want to get a single value from the data frame, we can provide an [index](https://swcarpentry.github.io/r-novice-inflammation/reference.html#index) in square brackets. The first number specifies the row and the second the column:

```r
dat[30, 20]
```

```r
[1] 16
```

Or we can indicate a continous rows and columns

```r
dat[1:4, 1:10]
```

```r
 V1 V2 V3 V4 V5 V6 V7 V8 V9 V10
1  0  0  1  3  1  2  4  7  8   3
2  0  1  2  1  2  1  3  2  2   6
3  0  1  1  3  3  2  6  2  5   9
4  0  0  2  0  4  2  2  1  6   7
```

With that in hand, let’s look at the help for `read.csv()`:

```r
?read.csv
```

There’s a lot of information there, but the most important part is the first couple of lines:

```r
read.csv(file, header = TRUE, sep = ",", quote = "\"",
         dec = ".", fill = TRUE, comment.char = "", ...)
```

This tells us that `read.csv()` has one argument, `file`, that doesn’t have a default value, and six others that do. At a minimum, we need to provide the file name. The are several other fucntions you can use to import such as `read.table()` and `read.delim()`. In fact, `read.csv()` and `read.delim()` are both wrappers with different arguments to conveniently import files with different delimiters. We. can view these changes by looking at the help menu:

```r
?read.delim
```

```r
read.delim(file, header = TRUE, sep = "\t", quote = "\"",
           dec = ".", fill = TRUE, comment.char = "", ...)
```

### Analyzing Multiple Datasets

We have created a function called `analyze` that creates graphs of the minimum, average, and maximum daily inflammation rates for a single data set:

```r
analyze <- function(filename) {
  # Plots the average, min, and max inflammation over time.
  # The argument or input is a character string representing the name and location of a CSV file.
  dat <- read.csv(file = filename, header = FALSE)
  avg_day_inflammation <- apply(dat, 2, mean)
  plot(avg_day_inflammation)
  max_day_inflammation <- apply(dat, 2, max)
  plot(max_day_inflammation)
  min_day_inflammation <- apply(dat, 2, min)
  plot(min_day_inflammation)
}

analyze("data/inflammation-01.csv")
```

Notice the structure of a function:

- **Function Name**: In this case the function is saved as `analyze`. We can call this function using `analyze()`
  
- **Formals**: This is the list of arguments that you can control when you call the function. In this case the only argument we can change when we call the function is the `filename`. The value of filename can be a path to a file that you want to analyze
  
- **Body**: This is the code within the brackets `{}`
  
- **Environment**: The data structure that finds the values assocaited each variable or names. Unless you specify a different environment, this will be the R Global Enviornment which is a list of objects created when we start an R Session.
  

**In the funtion above, what are we plotting? What functions are being called within our analyze function?What does the `2` in the nested function `apply` do?**

We can use it to analyze other data sets one by one:

```r
analyze("data/inflammation-02.csv")
```

### Processing multiple files using loops and lists

We can run the following command to get a list of all the .csv files in our directory.

```r
list.files(path = "data", pattern = "csv")
```

```r
 [1] "car-speeds-cleaned.csv" "car-speeds.csv"         "inflammation-01.csv"
 [4] "inflammation-02.csv"    "inflammation-03.csv"    "inflammation-04.csv"
 [7] "inflammation-05.csv"    "inflammation-06.csv"    "inflammation-07.csv"
[10] "inflammation-08.csv"    "inflammation-09.csv"    "inflammation-10.csv"
[13] "inflammation-11.csv"    "inflammation-12.csv"    "sample.csv"
[16] "small-01.csv"           "small-02.csv"           "small-03.csv"      
```

However, we only want the inflammation data files and none of the other text files. We can change the pattern to specify that we only want the inflammation data files.

```r
list.files(path = "data", pattern = "inflammation")
```

```r
 [1] "inflammation-01.csv"            "inflammation-02.csv"
 [3] "inflammation-03.csv"            "inflammation-04.csv"
 [5] "inflammation-05.csv"            "inflammation-06.csv"
 [7] "inflammation-07.csv"            "inflammation-08.csv"
 [9] "inflammation-09.csv"            "inflammation-10.csv"
[11] "inflammation-11.csv"            "inflammation-12.csv"
[13] "r-novice-inflammation-data.zip"
```

However, there is still a file present that we do not want `r-novice-inflammation-data.zip`. We can use expressions and wildcards to be more specific. Additionally, we can use the argument `full.names` to get the full relative path of our files.

```r
list.files(path = "data",  
          pattern = "inflammation-[0-9]{2}.csv",
          full.names = TRUE)
```

```r
 [1] "data/inflammation-01.csv"            "data/inflammation-02.csv"
 [3] "data/inflammation-03.csv"            "data/inflammation-04.csv"
 [5] "data/inflammation-05.csv"            "data/inflammation-06.csv"
 [7] "data/inflammation-07.csv"            "data/inflammation-08.csv"
 [9] "data/inflammation-09.csv"            "data/inflammation-10.csv"
[11] "data/inflammation-11.csv"            "data/inflammation-12.csv"
```

The pattern is composed of three parts:

- **Static**: This is the part of the filename that is always constant: `inflammation`
  
- **Variable**: This is the part of the filename that is different `[0-9]{2}` We are telling the function that our files have two digits and are numbered and range from 01 - 12.
  
- **Extension**: The `.csv` that tells us the specfic file types.
  

We can save this into a list called `filenames` and use that list to write a loop that will run our analyze function on each item listed.

```r
filenames <- list.files(path = "data",  
                        # Now follows a regular expression that matches:
                        pattern = "inflammation-[0-9]{2}.csv",
                        #          |            |        the standard file extension of comma-separated values
                        #          |            the variable parts (two digits, each between 0 and 9)
                        #          the static part of the filenames
                        full.names = TRUE)
filenames <- filenames[1:3]
for (f in filenames) {
  print(f)
  analyze(f)
}
```
