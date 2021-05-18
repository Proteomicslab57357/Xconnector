# Xconnector

Xconnector is a software package designed to easily retrieve, and visualize metabolomics data from different database sources. The goal of Xconnector is to connect different metabolomics databases in one place. The first five databases implemented in Xconnector are:

1. The Human Metabolome Database [HMDB](http://www.hmdb.ca/).
2. The Livestock Metabolome Database [LMDB](http://lmdb.ca/).
3. The Yeast Metabolome Database [YMDB](http://www.ymdb.ca/).
4. The Toxin and Toxin Target Database [T3DB](http://www.t3db.ca/).
5. ReSpect for Phytochemicals DataBase [ReSpect](http://spectra.psc.riken.jp/).
6. KEGG: Kyoto Encyclopedia of Genes and Genomes [KEGG](https://www.genome.jp/kegg/).
7. The Small Molecule Pathway Database [SMPDB](http://smpdb.ca/).

In future, we aim to include the most used databases for metabolites data.

# API & GUI Implementation

* The API function connects databases in Xconnector is made to be programmatically efficient. Using python generators implementation, only one query is called from the database each time by the API. This will reduce the memory used by Xconnector, as well as overcome the errors that could occur during the slow internet connection.

* After the API sends the output to the GUI. Xconnector utilises multithreading to allow efficient execution for the GUI, which allow multitasking and converting data between the GUI and the API.

# Workflow implementation

<img src="https://github.com/Proteomicslab57357/Xconnector/blob/master/image/Untitled%20Diagram%20(1).png" width="500" height="700">


* This workflow summarizes how Xconnector works. Using a graphical user interface (GUI), the user could import different type of input as database IDs or keyword from a CSV file or text file. Then send it to the application programming interface (API), the API connect one or all databases to retrieve and parse the information that made hits with the user input. Finally, the API send it back to the GUI to be displayed for the user.

# Software requirements

**A main requirement is a speed internet connection.**

* For Windows:

   For Windows 10 (64-bit system) users, two options to run Xconnector:
   
   1.1. Download it as normal Windows 10 software (executable program), with no need for any dependencies.
   
   1.2. Install it as a python package using pip. All the dependencies needed will be installed automatically (Note: Python >=3.7 is required).

* For Linux:
Install Xconnector as a python package using pip. All the dependencies needed will be installed automatically (Note: Python >=3.7 is required).

* For Mac:
Install Xconnector as a python package using pip. All the dependencies needed will be installed automatically (Note: Python >=3.7 is required).

# Software availability

* The Package can be downloaded using pip as: 

```python
pip install Xconnector
```
* For Windows 10 users the Executable program can be downloaded from: [Link](https://beta.57357.org/wp-content/themes/57357/programs/Xconnector.zip)

# How to download and install

* For Windows 10 users only.

  Download the zip file from [here](https://beta.57357.org/wp-content/themes/57357/programs/Xconnector.zip). After, extract it and double click on the icon named Xconnector.exe to run the software.

* For any operating system.

  Download and install Python >=3.7. Using pip download the Xconnector package as ( pip install Xconnector ). Then, Open python and run the package to run the GUI as:

```python
>>> from Xconnector import Xconnector
```
# Documentation

* For API function documentaion: Link

* For GUI documentation: [doc.](https://github.com/Proteomicslab57357/Xconnector/blob/master/Documentation.pdf)
