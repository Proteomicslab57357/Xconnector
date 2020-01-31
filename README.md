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
