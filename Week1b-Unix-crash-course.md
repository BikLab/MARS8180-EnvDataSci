# Crash course: UNIX Essentials

---

## Introduction to the UNIX Operating System

**By the end of this module you should be able to**

1. Identify the three major parts of the Unix Operating System
  
2. Explain how the Unix directory is structured
  
3. Demonstrate how to navigate the files and directory system using the Unix Shell
  

**Due to limited class time I have summarized the most important materials from these sources. On your own time, I highly recommend you read through the following tutorials**

- [Duke University Department of Computer Science](https://users.cs.duke.edu/~alvy/courses/unixtut/unixintro.html)
  
- [University of Georgia](https://wiki.gacrc.uga.edu/wiki/Training)
  
- [Software Carpentries](https://swcarpentry.github.io/shell-novice/index.html)
  

---

The UNIX Operating System was first developed in the 1960's and has been under constant development since then. There are many different versions of UNIX; however, the most popular varieties of UNIX are Linux and MacOS. Regardless of which version of UNIX you are using the operating system has the same three major components: the 1) Kernel, 2) Shell, 2) Programs.

## The Kernel

The kernel of UNIX is the hub of the operating system: it is the system that allocates time and memory to programs and handles the filestore and communications in response to system calls. 

## The Shell

The shell acts as an interface between the user and the kernel. When a user logs in, the login program checks the username and password, and then starts another program called the shell. The shell is a command line interpreter (CLI). It interprets the commands the user types in and arranges for them to be carried out. The commands are themselves programs: when they terminate, the shell gives the user another prompt (username@machinename> on our systems, but in this tutorial you will see % as the prompt in many examples...Do Not be Alarmed by this difference).

## Directory Structure

All files on your computer are organized heirarchly in a directory structure. At the top of the hierarchy is the "root" usually written as `/`. Personal directories are written in the `Users` folder. Many programs are written in `Bin` . The written locations are referred to as paths. The **absolute path** includes the entire path starting from the root directory while the **relative path** starts from a specific location. The paths can be thought of GPS directions with the. The aboslute path is the directions starting from your permanent residence, such as your house. But if need directions while you are bowling with your friend, you might require the relative path.

In the example below the aboslute path to the user *jsmith* directory is `/Users/jsmith/` Notice how the "root" is not written in the path - instead the front slash is used to indicate the root directory.

```mermaid
flowchart TD
    A[root] --> B[usr] 
    A[root] --> C[sys]
    A[root] --> D[var]
    A[root] --> E[temp] 
    A[root] --> F[Users]
    F[Users]--> G[desantiago]
    F[Users]--> H[Shared]
    G[jsmith] --> I[Desktop]
    G[jsmith] --> J[Documents]
    G[jsmith] --> K[Downloads]
    G[jsmith] --> L[Library]
    G[jsmith] --> M[Movies]
```

### Navigating files and directories

In order to orient yourself, you should first know where you are located in your computer. In order to **p**rint the **w**orking **d**irectory you can use the following command.

```bash
$ pwd
```

```bash
/Users/jsmith
```

You can **l**i**s**t the files and folders in your current working directory using the following command.

```bash
$ ls
```

```bash
Desktop    Downloads    Library    Pictures
Documents  Movies       Music      Public    
```

So far we know there are eight different folders or files in the directory. If we would like to classify the output we can use the flag `-F`

```bash
$ ls -F 
```

```bash
Desktop/    Downloads/    Library/    Pictures/
Documents/  Movies/       Music/      Public/    
```

The character after the output indicates what they are,

- a trailing `/` indicates that this is a directory
- `@` indicates a link
- `*` indicates an executable

In this case all the objects in the directory `/Users/jsmith/` are folders or subdirectories.

If you would like to navigate to a different directory you can use the **c**hange **d**irectories command

```bash
$ cd Desktop 
```

In this example, you only changed directories and DO NOT have an message outputed to your terminal. You can use the `pwd` command to print you new working directory.

```bash
$ pwd
```

```bash
/Users/jsmith/Desktop
```

If you would like to go back to your user directory you can either use the absolute path:

```bash
$ cd /Users/jsmith
```

or use a relative path:

```bash
$ cd ../
```

the `../` indicate the directory above the current directory.

### Accessing help menus or manuals

If you forget how to use a commnad there are several ways to pull up a help menu or the program's **man**ual.

```bash
$ man ls
```

### Syntax of a Shell Command

![](file:///Users/alejandrodesantiago/Downloads/shell_command_syntax.svg?msec=1735008803612)

The
